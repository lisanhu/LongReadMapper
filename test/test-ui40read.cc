//
// Created by lisanhu on 3/28/19.
//

#include <cstring>
#include "catch.hpp"
#include "../psascan/sa_use.h"




TEST_CASE("write ui40_t obj to disk", "[ui40_t]") {
    ui40_t num;
    num.high = 5;
    num.low = 0x01020304;
    FILE *fp = fopen("tmp", "w");
    fwrite(&num, sizeof(ui40_t), 1, fp);
    fclose(fp);
    fp = fopen("tmp", "r");
    auto buf = static_cast<ui40_t *>(malloc(sizeof(ui40_t) * 1));
    fread(buf, sizeof(ui40_t), 1, fp);
    fclose(fp);
    CHECK(buf->high == num.high);
    CHECK(buf->low == num.low);

    fp = fopen("tmp", "r");
    auto bytes = static_cast<char *>(malloc(sizeof(ui40_t) * 1));
    fread(bytes, 1, 8, fp);
    fclose(fp);
    for (int i = 0; i < 8; ++i) {
        printf("%d\n", bytes[i]);
    }

    ui40_t n = from_bytes((uint8_t *)bytes);
    CHECK(n.high == num.high);
    CHECK(n.low == num.low);
}

TEST_CASE("create .sa5 file with sa_use lib, and directly load ui40_t", "[ui40_t]") {
    const char *txt = "GAATTCCC";
    const size_t l = strlen(txt);
    FILE *fp = fopen("test.txt", "w");
    fwrite(txt, sizeof(char), l, fp);
    fclose(fp);
    sa_build("test.txt", 1 << 30);

    auto sas = static_cast<char *>(malloc(l * sizeof(ui40_t)));
    fp = fopen("test.txt.sa5", "r");
    size_t ll = fread(sas, 1, l * 5, fp);
    CHECK(ll == l * 5);
    fclose(fp);

    for (size_t i = 0; i < ll; ++i) {
//        printf("%d %d\n", sas[i].low, sas[i].high);
        printf("%d\n", sas[i]);
    }
}
