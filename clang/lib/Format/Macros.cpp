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

#include "clang/Basic/MemoryBufferCache.h"
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
        Diags(new DiagnosticIDs, new DiagnosticOptions, new IgnoringDiagConsumer),
        Headers(std::make_shared<HeaderSearchOptions>(), SourceMgr, Diags,
                LangOpts, nullptr),
        PP(makePPOpts(), Diags, LangOpts, SourceMgr, PCMCache, Headers,
           TheModuleLoader) {}

  std::shared_ptr<clang::PreprocessorOptions> makePPOpts() {
    auto PPOpts = std::make_shared<clang::PreprocessorOptions>();
    PPOpts->SingleFileParseMode = true;
    return PPOpts;
  }

  clang::SourceManager &SourceMgr;
  clang::LangOptions LangOpts;
  clang::DiagnosticsEngine Diags;
  clang::MemoryBufferCache PCMCache;
  clang::HeaderSearch Headers;
  clang::TrivialModuleLoader TheModuleLoader;
  clang::Preprocessor PP;
};

Macros::Macros(const std::vector<std::string> &Macros,
               clang::SourceManager &SourceMgr, const FormatStyle &Style)
    : PP(llvm::make_unique<PPState>(SourceMgr, Style)) {
  parseDefinitions(Macros);
}

Macros::~Macros() {}

void Macros::parseDefinitions(const std::vector<std::string> &Macros) {
  for (const std::string& Macro : Macros) {
    Definitions.push_back(llvm::MemoryBuffer::getMemBufferCopy("#define " + Macro + "\n", "<scratch space>"));
    process(Definitions.back().get());
  }
}

std::string Macros::process(llvm::MemoryBuffer *Buffer) {
  //clang::FileID FID = PP->SourceMgr.createFileID(std::move(Buffer));
  clang::FileID FID = PP->SourceMgr.createFileID(SourceManager::Unowned, Buffer);
  //clang::Lexer Lex(FID, Buffer.get(), PP->PP);
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
  return PP->PP.isMacroDefined(Name);
}

std::string Macros::Expand(llvm::StringRef Code) {
  auto Buffer = llvm::MemoryBuffer::getMemBufferCopy(Code, "<scratch space>");
  return process(Buffer.get());
}

}
}