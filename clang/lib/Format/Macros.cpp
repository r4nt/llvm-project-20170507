//===--- Macros.h - Format C++ code -----------------------------*- C++ -*-===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
///
/// \file
/// This file contains the implementation of Macros, which handles macro
/// configuration and expansion while formatting.
///
//===----------------------------------------------------------------------===//

#include "Macros.h"

#include "FormatToken.h"
#include "FormatTokenLexer.h"
#include "clang/Format/Format.h"
#include "clang/Lex/HeaderSearch.h"
#include "clang/Lex/HeaderSearchOptions.h"
#include "clang/Lex/Lexer.h"
#include "clang/Lex/ModuleLoader.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Lex/PreprocessorOptions.h"

namespace clang {
namespace format {

struct Macros::PPState {
  PPState(clang::SourceManager &SourceMgr, const FormatStyle &Style)
      : SourceMgr(SourceMgr), LangOpts(getFormattingLangOpts(Style)),
        Diags(new DiagnosticIDs, new DiagnosticOptions,
              new IgnoringDiagConsumer),
        Headers(std::make_shared<HeaderSearchOptions>(), SourceMgr, Diags,
                LangOpts, nullptr),
        PP(makePPOpts(), Diags, LangOpts, SourceMgr, Headers,
           TheModuleLoader) {}

  std::shared_ptr<clang::PreprocessorOptions> makePPOpts() {
    auto PPOpts = std::make_shared<clang::PreprocessorOptions>();
    PPOpts->SingleFileParseMode = true;
    return PPOpts;
  }

  clang::SourceManager &SourceMgr;
  clang::LangOptions LangOpts;
  clang::DiagnosticsEngine Diags;
  clang::HeaderSearch Headers;
  clang::TrivialModuleLoader TheModuleLoader;
  clang::Preprocessor PP;
};

struct Macros::Definition {
  StringRef Name;
  SmallVector<FormatToken *, 8> Params;
  SmallVector<FormatToken *, 8> Tokens;
};

class Macros::DefinitionParser {
public:
  DefinitionParser(ArrayRef<FormatToken *> Tokens) : Tokens(Tokens) {
    //for (const auto& T : Tokens) llvm::errs() << T->Tok.getName() << " ";
      //llvm::errs() << "\n";
    assert(!Tokens.empty());
    Tok = Tokens[0];
  }

  Macros::Definition parse() {
    if (!Tok->is(tok::identifier)) return {};
    Def.Name = Tok->TokenText;
    nextToken();
    if (Tok->is(tok::l_paren)) {
      parseParams();
    }  
    parseExpansion();
    //llvm::errs() << Def.Name << " / " << Def.Params.size() << " / " << Def.Tokens.size() << "\n";
    return Def;
  }

  void parseParams() {
    assert(Tok->is(tok::l_paren));
    nextToken();
    while (Tok->is(tok::identifier)) {
      Def.Params.push_back(Tok);
      //llvm::errs() << "PUSHING\n";
      nextToken();
      if (!Tok->is(tok::comma)) break;
      nextToken();
    }
    // FIXME: Error handling if not r_paren.
    nextToken();
  }

  void parseExpansion() {
    do {
      Def.Tokens.push_back(Tok);
      nextToken();
    } while(Tok->isNot(tok::eof));
    Def.Tokens.push_back(Tok);
  }  

  void nextToken() {
    if (I + 1 < Tokens.size())
      ++I;
    Tok = Tokens[I];
    //llvm::errs() << Tok->Tok.getName() << "\n";
    Tok->Finalized = true;
  }

private:
  size_t I = 0;
  FormatToken *Tok = nullptr;
  Definition Def;

  ArrayRef<FormatToken *> Tokens;
};

Macros::Macros(const std::vector<std::string> &Macros,
               clang::SourceManager &SourceMgr, const FormatStyle &Style,
               encoding::Encoding Encoding, 
               llvm::SpecificBumpPtrAllocator<FormatToken> &Allocator,
               IdentifierTable &IdentTable)
    : PP(llvm::make_unique<PPState>(SourceMgr, Style)), Style(Style),
      Encoding(Encoding), Allocator(Allocator), IdentTable(IdentTable) {
  parseDefinitions(Macros);
}

Macros::~Macros() {}

void Macros::parseDefinitions(const std::vector<std::string> &Macros) {
  for (const std::string& Macro : Macros) {
    /*Buffers.push_back(llvm::MemoryBuffer::getMemBufferCopy("#define " + Macro + "\n", "<scratch space>"));
    process(Buffers.back().get());
*/
    Buffers.push_back(llvm::MemoryBuffer::getMemBufferCopy(Macro, "<scratch space>"));
    clang::FileID FID = PP->SourceMgr.createFileID(SourceManager::Unowned, Buffers.back().get());
    FormatTokenLexer Lex(PP->SourceMgr, FID, 0, Style, Encoding, Allocator, IdentTable);
    DefinitionParser Parser(Lex.lex());
    auto Definition = Parser.parse();
    Definitions[Definition.Name] = Definition;
  }
}

std::string Macros::process(llvm::MemoryBuffer *Buffer) {
  //clang::FileID FID = PP->SourceMgr.createFileID(std::move(Buffer));
  //clang::Lexer Lex(FID, Buffer.get(), PP->PP);
  clang::FileID FID = PP->SourceMgr.createFileID(SourceManager::Unowned, Buffer);
  PP->PP.EnterSourceFile(FID, nullptr, SourceLocation());
  std::string Result;
  Token Tok;
  do {
    PP->PP.Lex(Tok);
    if (!Tok.is(tok::eof)) {
      if (!Result.empty()) Result.append(" ");
      Result.append(PP->SourceMgr.getCharacterData(Tok.getLocation()),
                    Tok.getLength());
    } 
  } while(!Tok.is(tok::eof));
  //llvm::errs() << "defined? " << PP->PP.isMacroDefined("CLASS")<< "\n";
  //llvm::errs() << "LETS SEE: " << Result << "\n";
  return Result;
}

bool Macros::Defined(llvm::StringRef Name) {
  return Definitions.find(Name) != Definitions.end();
  //return PP->PP.isMacroDefined(Name);
}
/*
std::string Macros::Expand(llvm::StringRef Code) {
  auto Buffer = llvm::MemoryBuffer::getMemBufferCopy(Code, "<scratch space>");
  return process(Buffer.get());
}
*/

llvm::SmallVector<FormatToken*, 8> Macros::Expand2(FormatToken *ID, llvm::ArrayRef<llvm::SmallVector<FormatToken *, 8>> Args) {
  /*llvm::errs() << "Name: " << Name << "\n";
  llvm::errs() << "Args: " << Args.size() << "\n";
  llvm::errs() << "Parm: " << Definitions[Name].Params.size() << "\n";
  llvm::errs() << "Expa: " << Definitions[Name].Tokens.size() << "\n";*/
  SmallVector<FormatToken*, 8> Result;
  const Definition &Def = Definitions[ID->TokenText];
  llvm::StringMap<int> ArgMap;
  for (int I = 0, E = Def.Params.size(); I != E; ++I) {
    ArgMap[Def.Params[I]->TokenText] = I;
  }
  bool First = true;
  for (FormatToken *Tok : Definitions[ID->TokenText].Tokens) {
    if (Tok->is(tok::identifier)) {
      auto I = ArgMap.find(Tok->TokenText);
      if (I != ArgMap.end()) {
        if (I->getValue() < Args.size()) {
          for (const auto &Tok : Args[I->getValue()]) {
            if (Tok->Macro == MS_None)
              Tok->Macro = MS_Expansion;
            Tok->ExpandedFrom.push_back(ID);
            if (First) {
              Tok->StartOfExpansion = true;
            }
            /*
            if (Tok->ExpandedFrom == nullptr) {
              Tok->ExpandedFrom = ID;
            }
            */
            Result.push_back(Tok);
            First = false;
          }
        }
        continue;
      }
    }
    FormatToken *New = new (Allocator.Allocate()) FormatToken;
    Tok->copyInto(*New);
    //New->ExpandedFrom = ID;
    assert(New->Macro == MS_None);
    New->Macro = MS_Hidden;
    New->ExpandedFrom.push_back(ID);
    if (First) {
      New->StartOfExpansion = true;
    }
    Result.push_back(New);
    First = false;
  }
  assert(Result.size() >= 2);
  ++Result[Result.size()-2]->EndOfExpansion;
  return Result;
}

//llvm::ArrayRef<FormatToken*> Macros::Expand(llvm::ArrayRef<FormatToken*> Call) {
  /*
  assert(!Call.empty());
  PP->PP.EnterSourceFile(PP->SourceMgr.getMainFileID(), nullptr, Call[0]->getStartOfNonWhitespace());
  llvm::SpecificBumpPtrAllocator<FormatToken> Allocator;
  FormatTokenLexer Lex(PP->SourceMgr, PP->SourceMgr.getMainFileID(),
                       0, Style, Encoding, Allocator);
  //ArrayRef<FormatToken*> Result = 
  //Lex.lex();
  auto End = [&](FormatToken *Tok) {
    llvm::errs() << "EEERRRRRRS\n";
    llvm::errs() << Tok->TokenText << "\n";
    return PP->SourceMgr.getExpansionLoc(Tok->getStartOfNonWhitespace()) !=
           PP->SourceMgr.getExpansionLoc(Call[0]->getStartOfNonWhitespace());
  };
  return Lex.lexUntil(End);
  */
//  return {};
//}

} // namespace format
} // namespace clang
