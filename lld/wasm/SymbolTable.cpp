//===- SymbolTable.cpp ----------------------------------------------------===//
//
//                             The LLVM Linker
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "SymbolTable.h"
#include "Config.h"
#include "InputChunks.h"
#include "InputGlobal.h"
#include "WriterUtils.h"
#include "lld/Common/ErrorHandler.h"
#include "lld/Common/Memory.h"
#include "llvm/ADT/SetVector.h"

#define DEBUG_TYPE "lld"

using namespace llvm;
using namespace llvm::wasm;
using namespace lld;
using namespace lld::wasm;

SymbolTable *lld::wasm::Symtab;

void SymbolTable::addFile(InputFile *File) {
  log("Processing: " + toString(File));
  File->parse();

  if (auto *F = dyn_cast<ObjFile>(File))
    ObjectFiles.push_back(F);
}

void SymbolTable::reportRemainingUndefines() {
  SetVector<Symbol *> Undefs;
  for (Symbol *Sym : SymVector) {
    if (Sym->isUndefined() && !Sym->isWeak() &&
        Config->AllowUndefinedSymbols.count(Sym->getName()) == 0) {
      Undefs.insert(Sym);
    }
  }

  if (Undefs.empty())
    return;

  for (ObjFile *File : ObjectFiles)
    for (Symbol *Sym : File->getSymbols())
      if (Undefs.count(Sym))
        error(toString(File) + ": undefined symbol: " + toString(*Sym));

  for (Symbol *Sym : Undefs)
    if (!Sym->getFile())
      error("undefined symbol: " + toString(*Sym));
}

Symbol *SymbolTable::find(StringRef Name) {
  auto It = SymMap.find(CachedHashStringRef(Name));
  if (It == SymMap.end())
    return nullptr;
  return It->second;
}

std::pair<Symbol *, bool> SymbolTable::insert(StringRef Name) {
  Symbol *&Sym = SymMap[CachedHashStringRef(Name)];
  if (Sym)
    return {Sym, false};
  Sym = reinterpret_cast<Symbol *>(make<SymbolUnion>());
  SymVector.emplace_back(Sym);
  return {Sym, true};
}

// Check the type of new symbol matches that of the symbol is replacing.
// For functions this can also involve verifying that the signatures match.
static void checkSymbolTypes(const Symbol &Existing, const InputFile &F,
                             WasmSymbolType NewType,
                             const WasmSignature *NewFunctionSig,
                             const WasmGlobalType *NewGlobalType) {
  if (Existing.isLazy())
    return;

  WasmSymbolType ExistingType = Existing.getWasmType();

  // First check the symbol types match (i.e. either both are function
  // symbols or both are data symbols).
  if (NewType != ExistingType) {
    error("symbol type mismatch: " + Existing.getName() + "\n>>> defined as " +
          toString(ExistingType) + " in " + toString(Existing.getFile()) +
          "\n>>> defined as " + toString(NewType) + " in " + F.getName());
    return;
  }

  // For function/global symbols, optionally check the type matches too.
  if (NewType == WASM_SYMBOL_TYPE_DATA || !Config->CheckSignatures)
    return;

  DEBUG(dbgs() << "checkSymbolTypes: " << Existing.getName() << "\n");

  auto ReportError = [&](const Twine &Old, const Twine &New) {
    error(toString(NewType) + " type mismatch: " + Existing.getName() +
          "\n>>> defined as " + Old + " in " + toString(Existing.getFile()) +
          "\n>>> defined as " + New + " in " + F.getName());
  };

  if (NewType == WASM_SYMBOL_TYPE_FUNCTION) {
    // Skip the signature check if the existing function has no signature (e.g.
    // if it is an undefined symbol generated by --undefined command line flag).
    auto &Sym = cast<FunctionSymbol>(Existing);
    const WasmSignature *OldSig = Sym.getFunctionType();
    if (!OldSig)
      return;

    assert(NewFunctionSig);
    if (*NewFunctionSig == *OldSig)
      return;

    ReportError(toString(*OldSig), toString(*NewFunctionSig));
  } else {
    auto &Sym = cast<GlobalSymbol>(Existing);

    assert(NewGlobalType != nullptr);
    const WasmGlobalType *OldType = Sym.getGlobalType();
    if (*NewGlobalType == *OldType)
      return;

    ReportError(toString(*OldType), toString(*NewGlobalType));
  }
}

DefinedFunction *SymbolTable::addSyntheticFunction(StringRef Name,
                                                   const WasmSignature *Type,
                                                   uint32_t Flags) {
  DEBUG(dbgs() << "addSyntheticFunction: " << Name << "\n");
  Symbol *S;
  bool WasInserted;
  std::tie(S, WasInserted) = insert(Name);
  assert(WasInserted);
  return replaceSymbol<DefinedFunction>(S, Name, Flags, Type);
}

DefinedData *SymbolTable::addSyntheticDataSymbol(StringRef Name,
                                                 uint32_t Flags) {
  DEBUG(dbgs() << "addSyntheticDataSymbol: " << Name << "\n");
  Symbol *S;
  bool WasInserted;
  std::tie(S, WasInserted) = insert(Name);
  assert(WasInserted);
  return replaceSymbol<DefinedData>(S, Name, Flags);
}

DefinedGlobal *SymbolTable::addSyntheticGlobal(StringRef Name, uint32_t Flags,
                                               InputGlobal *Global) {
  DEBUG(dbgs() << "addSyntheticGlobal: " << Name << " -> " << Global << "\n");
  Symbol *S;
  bool WasInserted;
  std::tie(S, WasInserted) = insert(Name);
  assert(WasInserted);
  return replaceSymbol<DefinedGlobal>(S, Name, Flags, nullptr, Global);
}

static bool shouldReplace(const Symbol &Existing, InputFile *NewFile,
                          WasmSymbolType NewType, uint32_t NewFlags,
                          const WasmSignature *NewFuncType = nullptr,
                          const WasmGlobalType *NewGlobalType = nullptr) {

  // If existing symbol is lazy, replace it without checking types since
  // lazy symbols don't have any type information.
  if (Existing.isLazy()) {
    DEBUG(dbgs() << "replacing existing lazy symbol: " << Existing.getName()
                 << "\n");
    return true;
  }

  // Now we have two wasm symbols, and all wasm symbols that have the same
  // symbol name must have the same type, even if they are undefined. This
  // is different from ELF because symbol types are not that significant
  // in ELF, and undefined symbols in ELF don't have type in the first place.
  checkSymbolTypes(Existing, *NewFile, NewType, NewFuncType, NewGlobalType);

  // If existing symbol is undefined, replace it.
  if (!Existing.isDefined()) {
    DEBUG(dbgs() << "resolving existing undefined symbol: "
                 << Existing.getName() << "\n");
    return true;
  }

  // Now we have two defined symbols. If the new one is weak, we can ignore it.
  if ((NewFlags & WASM_SYMBOL_BINDING_MASK) == WASM_SYMBOL_BINDING_WEAK) {
    DEBUG(dbgs() << "existing symbol takes precedence\n");
    return false;
  }

  // If the existing symbol is weak, we should replace it.
  if (Existing.isWeak()) {
    DEBUG(dbgs() << "replacing existing weak symbol\n");
    return true;
  }

  // Neither symbol is week. They conflict.
  error("duplicate symbol: " + toString(Existing) + "\n>>> defined in " +
        toString(Existing.getFile()) + "\n>>> defined in " + toString(NewFile));
  return true;
}

Symbol *SymbolTable::addDefinedFunction(StringRef Name, uint32_t Flags,
                                        InputFile *F, InputFunction *Function) {
  DEBUG(dbgs() << "addDefinedFunction: " << Name << "\n");
  Symbol *S;
  bool WasInserted;
  std::tie(S, WasInserted) = insert(Name);
  if (WasInserted || shouldReplace(*S, F, WASM_SYMBOL_TYPE_FUNCTION, Flags,
                                   &Function->Signature))
    replaceSymbol<DefinedFunction>(S, Name, Flags, F, Function);
  return S;
}

Symbol *SymbolTable::addDefinedData(StringRef Name, uint32_t Flags,
                                    InputFile *F, InputSegment *Segment,
                                    uint32_t Address, uint32_t Size) {
  DEBUG(dbgs() << "addDefinedData:" << Name << " addr:" << Address << "\n");
  Symbol *S;
  bool WasInserted;
  std::tie(S, WasInserted) = insert(Name);
  if (WasInserted || shouldReplace(*S, F, WASM_SYMBOL_TYPE_DATA, Flags))
    replaceSymbol<DefinedData>(S, Name, Flags, F, Segment, Address, Size);
  return S;
}

Symbol *SymbolTable::addDefinedGlobal(StringRef Name, uint32_t Flags,
                                      InputFile *F, InputGlobal *Global) {
  DEBUG(dbgs() << "addDefinedGlobal:" << Name << "\n");
  Symbol *S;
  bool WasInserted;
  std::tie(S, WasInserted) = insert(Name);
  if (WasInserted || shouldReplace(*S, F, WASM_SYMBOL_TYPE_GLOBAL, Flags,
                                   nullptr, &Global->getType()))
    replaceSymbol<DefinedGlobal>(S, Name, Flags, F, Global);
  return S;
}

Symbol *SymbolTable::addUndefined(StringRef Name, WasmSymbolType Type,
                                  uint32_t Flags, InputFile *F,
                                  const WasmSignature *FunctionType,
                                  const WasmGlobalType *GlobalType) {
  DEBUG(dbgs() << "addUndefined type=" << Type << ": " << Name << "\n");

  Symbol *S;
  bool WasInserted;
  std::tie(S, WasInserted) = insert(Name);

  if (WasInserted) {
    switch (Type) {
    case WASM_SYMBOL_TYPE_FUNCTION:
      replaceSymbol<UndefinedFunction>(S, Name, Flags, F, FunctionType);
      break;
    case WASM_SYMBOL_TYPE_GLOBAL:
      replaceSymbol<UndefinedGlobal>(S, Name, Flags, F, GlobalType);
      break;
    case WASM_SYMBOL_TYPE_DATA:
      replaceSymbol<UndefinedData>(S, Name, Flags, F);
      break;
    }
    return S;
  }

  if (auto *Lazy = dyn_cast<LazySymbol>(S)) {
    DEBUG(dbgs() << "resolved by existing lazy\n");
    cast<ArchiveFile>(Lazy->getFile())->addMember(&Lazy->getArchiveSymbol());
    return S;
  }

  if (S->isDefined()) {
    DEBUG(dbgs() << "resolved by existing\n");
    checkSymbolTypes(*S, *F, Type, FunctionType, GlobalType);
  }

  return S;
}

void SymbolTable::addLazy(ArchiveFile *F, const Archive::Symbol *Sym) {
  DEBUG(dbgs() << "addLazy: " << Sym->getName() << "\n");
  StringRef Name = Sym->getName();

  Symbol *S;
  bool WasInserted;
  std::tie(S, WasInserted) = insert(Name);

  if (WasInserted) {
    replaceSymbol<LazySymbol>(S, Name, F, *Sym);
    return;
  }

  // If there is an existing undefined symbol, load a new one from the archive.
  if (S->isUndefined()) {
    DEBUG(dbgs() << "replacing existing undefined\n");
    F->addMember(Sym);
  }
}

bool SymbolTable::addComdat(StringRef Name, ObjFile *F) {
  DEBUG(dbgs() << "addComdat: " << Name << "\n");
  ObjFile *&File = ComdatMap[CachedHashStringRef(Name)];
  if (File) {
    DEBUG(dbgs() << "COMDAT already defined\n");
    return false;
  }
  File = F;
  return true;
}

ObjFile *SymbolTable::findComdat(StringRef Name) const {
  auto It = ComdatMap.find(CachedHashStringRef(Name));
  return It == ComdatMap.end() ? nullptr : It->second;
}
