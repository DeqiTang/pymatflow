#include <gtest/gtest.h>

#include <atomsciflow/utils.h>

TEST(utils, version) {
    EXPECT_EQ(atomsciflow::version(), "0.0.0");
}
