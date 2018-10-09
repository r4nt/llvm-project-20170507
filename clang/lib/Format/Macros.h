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

#include "Encoding.h"
#include "llvm/ADT/StringRef.h"
#include "llvm/ADT/ArrayRef.h"

namespace llvm { class MemoryBuffer; }

namespace clang {
class IdentifierTable;
class SourceManager;

namespace format {
struct FormatStyle;
struct FormatToken;

class Macros {
public:
  Macros(const std::vector<std::string> &Macros,
         clang::SourceManager &SourceMgr, const FormatStyle &Style,
         encoding::Encoding encoding,
         llvm::SpecificBumpPtrAllocator<FormatToken> &Allocator,
         IdentifierTable &IdentTable);
  ~Macros();
  bool Defined(llvm::StringRef Name);
  //std::string Expand(llvm::StringRef Code);

  llvm::SmallVector<FormatToken *, 8>
  Expand2(FormatToken *ID,
         llvm::ArrayRef<llvm::SmallVector<FormatToken *, 8>> Args);
  //llvm::ArrayRef<FormatToken*> Expand(llvm::ArrayRef<FormatToken*> Call);

private:
  struct PPState;
  struct Definition;
  class  DefinitionParser;

  void parseDefinitions(const std::vector<std::string> &Macros);
  std::string process(llvm::MemoryBuffer *Buffer);

  std::unique_ptr<PPState> PP;
  const FormatStyle &Style;
  encoding::Encoding Encoding;
  llvm::SpecificBumpPtrAllocator<FormatToken> &Allocator;
  IdentifierTable &IdentTable;
  std::vector<std::unique_ptr<llvm::MemoryBuffer>> Buffers;

  llvm::StringMap<Definition> Definitions;
};

} // namespace format
} // namespace clang

#endif
