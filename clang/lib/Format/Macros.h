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
/// This file contains the declaration of Macros, which handles macro
/// configuration and expansion while formatting.
///
//===----------------------------------------------------------------------===//

#ifndef LLVM_CLANG_LIB_FORMAT_MACROS_H
#define LLVM_CLANG_LIB_FORMAT_MACROS_H

#include <string>
#include <vector>
#include <unordered_map>

#include "llvm/ADT/StringRef.h"

namespace llvm { class MemoryBuffer; }

namespace clang {
class SourceManager;

namespace format {
struct FormatStyle;

class Macros {
public:
  Macros(const std::vector<std::string> &Macros,
         clang::SourceManager &SourceMgr, const FormatStyle &Style);
  ~Macros();
  bool Defined(llvm::StringRef Name);
  std::string Expand(llvm::StringRef Code);

private:
  struct PPState;

  void parseDefinitions(const std::vector<std::string> &Macros);
  std::string process(llvm::MemoryBuffer *Buffer);

  std::unique_ptr<PPState> PP;
  std::vector<std::unique_ptr<llvm::MemoryBuffer>> Definitions;
};

} // namespace format
} // namespace clang

#endif
