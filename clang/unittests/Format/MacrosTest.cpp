#include "../../lib/Format/Macros.h"

#include "clang/Basic/SourceManager.h"
#include "clang/Format/Format.h"

#include "gtest/gtest.h"

namespace clang {
namespace format {
namespace {

/*
TEST(MacrosTest, ReplacesIdentifier) {
  SourceManagerForFile SourceMgr("test", "");
  Macros M(std::vector<std::string>({"X x"}), SourceMgr.get(), getLLVMStyle());
  EXPECT_TRUE(M.Defined("X"));
  EXPECT_EQ("x", M.Expand("X"));
}

TEST(MacrosTest, ReplacesCall) {
  SourceManagerForFile SourceMgr("test", "");
  Macros M(std::vector<std::string>({"ID(x) x"}), SourceMgr.get(), getLLVMStyle());
  EXPECT_TRUE(M.Defined("ID"));
  EXPECT_EQ("a * b", M.Expand("ID(a * b)"));
}
*/

} // namespace
} // namespace format
} // namespace clang
