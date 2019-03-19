; NOTE: Assertions have been autogenerated by utils/update_llc_test_checks.py
; RUN: llc < %s -mtriple=powerpc64le-- -verify-machineinstrs | FileCheck %s

; There are at least 3 potential patterns corresponding to an unsigned saturated add: min, cmp with sum, cmp with not.
; Test each of those patterns with i8/i16/i32/i64.
; Test each of those with a constant operand and a variable operand.
; Test each of those with a 128-bit vector type.

define i8 @unsigned_sat_constant_i8_using_min(i8 %x) {
; CHECK-LABEL: unsigned_sat_constant_i8_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    clrlwi 5, 3, 24
; CHECK-NEXT:    li 4, -43
; CHECK-NEXT:    cmplwi 5, 213
; CHECK-NEXT:    isel 3, 3, 4, 0
; CHECK-NEXT:    addi 3, 3, 42
; CHECK-NEXT:    blr
  %c = icmp ult i8 %x, -43
  %s = select i1 %c, i8 %x, i8 -43
  %r = add i8 %s, 42
  ret i8 %r
}

define i8 @unsigned_sat_constant_i8_using_cmp_sum(i8 %x) {
; CHECK-LABEL: unsigned_sat_constant_i8_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    rlwinm 3, 3, 0, 24, 31
; CHECK-NEXT:    addi 3, 3, 42
; CHECK-NEXT:    andi. 4, 3, 256
; CHECK-NEXT:    li 4, -1
; CHECK-NEXT:    isel 3, 3, 4, 2
; CHECK-NEXT:    blr
  %a = add i8 %x, 42
  %c = icmp ugt i8 %x, %a
  %r = select i1 %c, i8 -1, i8 %a
  ret i8 %r
}

define i8 @unsigned_sat_constant_i8_using_cmp_notval(i8 %x) {
; CHECK-LABEL: unsigned_sat_constant_i8_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    clrlwi 5, 3, 24
; CHECK-NEXT:    li 4, -1
; CHECK-NEXT:    addi 3, 3, 42
; CHECK-NEXT:    cmplwi 5, 213
; CHECK-NEXT:    isel 3, 4, 3, 1
; CHECK-NEXT:    blr
  %a = add i8 %x, 42
  %c = icmp ugt i8 %x, -43
  %r = select i1 %c, i8 -1, i8 %a
  ret i8 %r
}

define i16 @unsigned_sat_constant_i16_using_min(i16 %x) {
; CHECK-LABEL: unsigned_sat_constant_i16_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    clrlwi 5, 3, 16
; CHECK-NEXT:    li 4, -43
; CHECK-NEXT:    cmplwi 5, 65493
; CHECK-NEXT:    isel 3, 3, 4, 0
; CHECK-NEXT:    addi 3, 3, 42
; CHECK-NEXT:    blr
  %c = icmp ult i16 %x, -43
  %s = select i1 %c, i16 %x, i16 -43
  %r = add i16 %s, 42
  ret i16 %r
}

define i16 @unsigned_sat_constant_i16_using_cmp_sum(i16 %x) {
; CHECK-LABEL: unsigned_sat_constant_i16_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    rlwinm 3, 3, 0, 16, 31
; CHECK-NEXT:    addi 3, 3, 42
; CHECK-NEXT:    andis. 4, 3, 1
; CHECK-NEXT:    li 4, -1
; CHECK-NEXT:    isel 3, 3, 4, 2
; CHECK-NEXT:    blr
  %a = add i16 %x, 42
  %c = icmp ugt i16 %x, %a
  %r = select i1 %c, i16 -1, i16 %a
  ret i16 %r
}

define i16 @unsigned_sat_constant_i16_using_cmp_notval(i16 %x) {
; CHECK-LABEL: unsigned_sat_constant_i16_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    clrlwi 5, 3, 16
; CHECK-NEXT:    li 4, -1
; CHECK-NEXT:    addi 3, 3, 42
; CHECK-NEXT:    cmplwi 5, 65493
; CHECK-NEXT:    isel 3, 4, 3, 1
; CHECK-NEXT:    blr
  %a = add i16 %x, 42
  %c = icmp ugt i16 %x, -43
  %r = select i1 %c, i16 -1, i16 %a
  ret i16 %r
}

define i32 @unsigned_sat_constant_i32_using_min(i32 %x) {
; CHECK-LABEL: unsigned_sat_constant_i32_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    li 4, -43
; CHECK-NEXT:    cmplw 3, 4
; CHECK-NEXT:    isel 3, 3, 4, 0
; CHECK-NEXT:    addi 3, 3, 42
; CHECK-NEXT:    blr
  %c = icmp ult i32 %x, -43
  %s = select i1 %c, i32 %x, i32 -43
  %r = add i32 %s, 42
  ret i32 %r
}

define i32 @unsigned_sat_constant_i32_using_cmp_sum(i32 %x) {
; CHECK-LABEL: unsigned_sat_constant_i32_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addi 5, 3, 42
; CHECK-NEXT:    li 4, -1
; CHECK-NEXT:    cmplw 0, 5, 3
; CHECK-NEXT:    isel 3, 4, 5, 0
; CHECK-NEXT:    blr
  %a = add i32 %x, 42
  %c = icmp ugt i32 %x, %a
  %r = select i1 %c, i32 -1, i32 %a
  ret i32 %r
}

define i32 @unsigned_sat_constant_i32_using_cmp_notval(i32 %x) {
; CHECK-LABEL: unsigned_sat_constant_i32_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    li 4, -43
; CHECK-NEXT:    addi 5, 3, 42
; CHECK-NEXT:    cmplw 0, 3, 4
; CHECK-NEXT:    li 3, -1
; CHECK-NEXT:    isel 3, 3, 5, 1
; CHECK-NEXT:    blr
  %a = add i32 %x, 42
  %c = icmp ugt i32 %x, -43
  %r = select i1 %c, i32 -1, i32 %a
  ret i32 %r
}

define i64 @unsigned_sat_constant_i64_using_min(i64 %x) {
; CHECK-LABEL: unsigned_sat_constant_i64_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    li 4, -43
; CHECK-NEXT:    cmpld 3, 4
; CHECK-NEXT:    isel 3, 3, 4, 0
; CHECK-NEXT:    addi 3, 3, 42
; CHECK-NEXT:    blr
  %c = icmp ult i64 %x, -43
  %s = select i1 %c, i64 %x, i64 -43
  %r = add i64 %s, 42
  ret i64 %r
}

define i64 @unsigned_sat_constant_i64_using_cmp_sum(i64 %x) {
; CHECK-LABEL: unsigned_sat_constant_i64_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addi 5, 3, 42
; CHECK-NEXT:    li 4, -1
; CHECK-NEXT:    cmpld 5, 3
; CHECK-NEXT:    isel 3, 4, 5, 0
; CHECK-NEXT:    blr
  %a = add i64 %x, 42
  %c = icmp ugt i64 %x, %a
  %r = select i1 %c, i64 -1, i64 %a
  ret i64 %r
}

define i64 @unsigned_sat_constant_i64_using_cmp_notval(i64 %x) {
; CHECK-LABEL: unsigned_sat_constant_i64_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    li 4, -43
; CHECK-NEXT:    addi 5, 3, 42
; CHECK-NEXT:    cmpld 3, 4
; CHECK-NEXT:    li 3, -1
; CHECK-NEXT:    isel 3, 3, 5, 1
; CHECK-NEXT:    blr
  %a = add i64 %x, 42
  %c = icmp ugt i64 %x, -43
  %r = select i1 %c, i64 -1, i64 %a
  ret i64 %r
}

define i8 @unsigned_sat_variable_i8_using_min(i8 %x, i8 %y) {
; CHECK-LABEL: unsigned_sat_variable_i8_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    nor 5, 4, 4
; CHECK-NEXT:    clrlwi 6, 3, 24
; CHECK-NEXT:    clrlwi 7, 5, 24
; CHECK-NEXT:    cmplw 6, 7
; CHECK-NEXT:    isel 3, 3, 5, 0
; CHECK-NEXT:    add 3, 3, 4
; CHECK-NEXT:    blr
  %noty = xor i8 %y, -1
  %c = icmp ult i8 %x, %noty
  %s = select i1 %c, i8 %x, i8 %noty
  %r = add i8 %s, %y
  ret i8 %r
}

define i8 @unsigned_sat_variable_i8_using_cmp_sum(i8 %x, i8 %y) {
; CHECK-LABEL: unsigned_sat_variable_i8_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    rlwinm 4, 4, 0, 24, 31
; CHECK-NEXT:    rlwinm 3, 3, 0, 24, 31
; CHECK-NEXT:    add 3, 3, 4
; CHECK-NEXT:    andi. 4, 3, 256
; CHECK-NEXT:    li 4, -1
; CHECK-NEXT:    isel 3, 3, 4, 2
; CHECK-NEXT:    blr
  %a = add i8 %x, %y
  %c = icmp ugt i8 %x, %a
  %r = select i1 %c, i8 -1, i8 %a
  ret i8 %r
}

define i8 @unsigned_sat_variable_i8_using_cmp_notval(i8 %x, i8 %y) {
; CHECK-LABEL: unsigned_sat_variable_i8_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    nor 6, 4, 4
; CHECK-NEXT:    clrlwi 7, 3, 24
; CHECK-NEXT:    li 5, -1
; CHECK-NEXT:    add 3, 3, 4
; CHECK-NEXT:    clrlwi 6, 6, 24
; CHECK-NEXT:    cmplw 7, 6
; CHECK-NEXT:    isel 3, 5, 3, 1
; CHECK-NEXT:    blr
  %noty = xor i8 %y, -1
  %a = add i8 %x, %y
  %c = icmp ugt i8 %x, %noty
  %r = select i1 %c, i8 -1, i8 %a
  ret i8 %r
}

define i16 @unsigned_sat_variable_i16_using_min(i16 %x, i16 %y) {
; CHECK-LABEL: unsigned_sat_variable_i16_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    nor 5, 4, 4
; CHECK-NEXT:    clrlwi 6, 3, 16
; CHECK-NEXT:    clrlwi 7, 5, 16
; CHECK-NEXT:    cmplw 6, 7
; CHECK-NEXT:    isel 3, 3, 5, 0
; CHECK-NEXT:    add 3, 3, 4
; CHECK-NEXT:    blr
  %noty = xor i16 %y, -1
  %c = icmp ult i16 %x, %noty
  %s = select i1 %c, i16 %x, i16 %noty
  %r = add i16 %s, %y
  ret i16 %r
}

define i16 @unsigned_sat_variable_i16_using_cmp_sum(i16 %x, i16 %y) {
; CHECK-LABEL: unsigned_sat_variable_i16_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    rlwinm 4, 4, 0, 16, 31
; CHECK-NEXT:    rlwinm 3, 3, 0, 16, 31
; CHECK-NEXT:    add 3, 3, 4
; CHECK-NEXT:    andis. 4, 3, 1
; CHECK-NEXT:    li 4, -1
; CHECK-NEXT:    isel 3, 3, 4, 2
; CHECK-NEXT:    blr
  %a = add i16 %x, %y
  %c = icmp ugt i16 %x, %a
  %r = select i1 %c, i16 -1, i16 %a
  ret i16 %r
}

define i16 @unsigned_sat_variable_i16_using_cmp_notval(i16 %x, i16 %y) {
; CHECK-LABEL: unsigned_sat_variable_i16_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    nor 6, 4, 4
; CHECK-NEXT:    clrlwi 7, 3, 16
; CHECK-NEXT:    li 5, -1
; CHECK-NEXT:    add 3, 3, 4
; CHECK-NEXT:    clrlwi 6, 6, 16
; CHECK-NEXT:    cmplw 7, 6
; CHECK-NEXT:    isel 3, 5, 3, 1
; CHECK-NEXT:    blr
  %noty = xor i16 %y, -1
  %a = add i16 %x, %y
  %c = icmp ugt i16 %x, %noty
  %r = select i1 %c, i16 -1, i16 %a
  ret i16 %r
}

define i32 @unsigned_sat_variable_i32_using_min(i32 %x, i32 %y) {
; CHECK-LABEL: unsigned_sat_variable_i32_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    nor 5, 4, 4
; CHECK-NEXT:    cmplw 3, 5
; CHECK-NEXT:    isel 3, 3, 5, 0
; CHECK-NEXT:    add 3, 3, 4
; CHECK-NEXT:    blr
  %noty = xor i32 %y, -1
  %c = icmp ult i32 %x, %noty
  %s = select i1 %c, i32 %x, i32 %noty
  %r = add i32 %s, %y
  ret i32 %r
}

define i32 @unsigned_sat_variable_i32_using_cmp_sum(i32 %x, i32 %y) {
; CHECK-LABEL: unsigned_sat_variable_i32_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    add 4, 3, 4
; CHECK-NEXT:    li 5, -1
; CHECK-NEXT:    cmplw 0, 4, 3
; CHECK-NEXT:    isel 3, 5, 4, 0
; CHECK-NEXT:    blr
  %a = add i32 %x, %y
  %c = icmp ugt i32 %x, %a
  %r = select i1 %c, i32 -1, i32 %a
  ret i32 %r
}

define i32 @unsigned_sat_variable_i32_using_cmp_notval(i32 %x, i32 %y) {
; CHECK-LABEL: unsigned_sat_variable_i32_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    nor 6, 4, 4
; CHECK-NEXT:    li 5, -1
; CHECK-NEXT:    cmplw 3, 6
; CHECK-NEXT:    add 3, 3, 4
; CHECK-NEXT:    isel 3, 5, 3, 1
; CHECK-NEXT:    blr
  %noty = xor i32 %y, -1
  %a = add i32 %x, %y
  %c = icmp ugt i32 %x, %noty
  %r = select i1 %c, i32 -1, i32 %a
  ret i32 %r
}

define i64 @unsigned_sat_variable_i64_using_min(i64 %x, i64 %y) {
; CHECK-LABEL: unsigned_sat_variable_i64_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    not 5, 4
; CHECK-NEXT:    cmpld 3, 5
; CHECK-NEXT:    isel 3, 3, 5, 0
; CHECK-NEXT:    add 3, 3, 4
; CHECK-NEXT:    blr
  %noty = xor i64 %y, -1
  %c = icmp ult i64 %x, %noty
  %s = select i1 %c, i64 %x, i64 %noty
  %r = add i64 %s, %y
  ret i64 %r
}

define i64 @unsigned_sat_variable_i64_using_cmp_sum(i64 %x, i64 %y) {
; CHECK-LABEL: unsigned_sat_variable_i64_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    add 4, 3, 4
; CHECK-NEXT:    li 5, -1
; CHECK-NEXT:    cmpld 4, 3
; CHECK-NEXT:    isel 3, 5, 4, 0
; CHECK-NEXT:    blr
  %a = add i64 %x, %y
  %c = icmp ugt i64 %x, %a
  %r = select i1 %c, i64 -1, i64 %a
  ret i64 %r
}

define i64 @unsigned_sat_variable_i64_using_cmp_notval(i64 %x, i64 %y) {
; CHECK-LABEL: unsigned_sat_variable_i64_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    not 6, 4
; CHECK-NEXT:    li 5, -1
; CHECK-NEXT:    cmpld 3, 6
; CHECK-NEXT:    add 3, 3, 4
; CHECK-NEXT:    isel 3, 5, 3, 1
; CHECK-NEXT:    blr
  %noty = xor i64 %y, -1
  %a = add i64 %x, %y
  %c = icmp ugt i64 %x, %noty
  %r = select i1 %c, i64 -1, i64 %a
  ret i64 %r
}

define <16 x i8> @unsigned_sat_constant_v16i8_using_min(<16 x i8> %x) {
; CHECK-LABEL: unsigned_sat_constant_v16i8_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI24_0@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI24_0@toc@l
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    addis 3, 2, .LCPI24_1@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI24_1@toc@l
; CHECK-NEXT:    vcmpgtub 4, 3, 2
; CHECK-NEXT:    xxsel 34, 35, 34, 36
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    vaddubm 2, 2, 3
; CHECK-NEXT:    blr
  %c = icmp ult <16 x i8> %x, <i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43>
  %s = select <16 x i1> %c, <16 x i8> %x, <16 x i8> <i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43>
  %r = add <16 x i8> %s, <i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42>
  ret <16 x i8> %r
}

define <16 x i8> @unsigned_sat_constant_v16i8_using_cmp_sum(<16 x i8> %x) {
; CHECK-LABEL: unsigned_sat_constant_v16i8_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI25_0@toc@ha
; CHECK-NEXT:    vspltisb 4, -1
; CHECK-NEXT:    addi 3, 3, .LCPI25_0@toc@l
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    vaddubm 3, 2, 3
; CHECK-NEXT:    vcmpgtub 2, 2, 3
; CHECK-NEXT:    xxsel 34, 35, 36, 34
; CHECK-NEXT:    blr
  %a = add <16 x i8> %x, <i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42>
  %c = icmp ugt <16 x i8> %x, %a
  %r = select <16 x i1> %c, <16 x i8> <i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1>, <16 x i8> %a
  ret <16 x i8> %r
}

define <16 x i8> @unsigned_sat_constant_v16i8_using_cmp_notval(<16 x i8> %x) {
; CHECK-LABEL: unsigned_sat_constant_v16i8_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI26_1@toc@ha
; CHECK-NEXT:    vspltisb 5, -1
; CHECK-NEXT:    addi 3, 3, .LCPI26_1@toc@l
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    addis 3, 2, .LCPI26_0@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI26_0@toc@l
; CHECK-NEXT:    vcmpgtub 3, 2, 3
; CHECK-NEXT:    lvx 4, 0, 3
; CHECK-NEXT:    vaddubm 2, 2, 4
; CHECK-NEXT:    xxsel 34, 34, 37, 35
; CHECK-NEXT:    blr
  %a = add <16 x i8> %x, <i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42, i8 42>
  %c = icmp ugt <16 x i8> %x, <i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43, i8 -43>
  %r = select <16 x i1> %c, <16 x i8> <i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1>, <16 x i8> %a
  ret <16 x i8> %r
}

define <8 x i16> @unsigned_sat_constant_v8i16_using_min(<8 x i16> %x) {
; CHECK-LABEL: unsigned_sat_constant_v8i16_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI27_0@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI27_0@toc@l
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    addis 3, 2, .LCPI27_1@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI27_1@toc@l
; CHECK-NEXT:    vcmpgtuh 4, 3, 2
; CHECK-NEXT:    xxsel 34, 35, 34, 36
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    vadduhm 2, 2, 3
; CHECK-NEXT:    blr
  %c = icmp ult <8 x i16> %x, <i16 -43, i16 -43, i16 -43, i16 -43, i16 -43, i16 -43, i16 -43, i16 -43>
  %s = select <8 x i1> %c, <8 x i16> %x, <8 x i16> <i16 -43, i16 -43, i16 -43, i16 -43, i16 -43, i16 -43, i16 -43, i16 -43>
  %r = add <8 x i16> %s, <i16 42, i16 42, i16 42, i16 42, i16 42, i16 42, i16 42, i16 42>
  ret <8 x i16> %r
}

define <8 x i16> @unsigned_sat_constant_v8i16_using_cmp_sum(<8 x i16> %x) {
; CHECK-LABEL: unsigned_sat_constant_v8i16_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI28_0@toc@ha
; CHECK-NEXT:    vspltisb 4, -1
; CHECK-NEXT:    addi 3, 3, .LCPI28_0@toc@l
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    vadduhm 3, 2, 3
; CHECK-NEXT:    vcmpgtuh 2, 2, 3
; CHECK-NEXT:    xxsel 34, 35, 36, 34
; CHECK-NEXT:    blr
  %a = add <8 x i16> %x, <i16 42, i16 42, i16 42, i16 42, i16 42, i16 42, i16 42, i16 42>
  %c = icmp ugt <8 x i16> %x, %a
  %r = select <8 x i1> %c, <8 x i16> <i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1>, <8 x i16> %a
  ret <8 x i16> %r
}

define <8 x i16> @unsigned_sat_constant_v8i16_using_cmp_notval(<8 x i16> %x) {
; CHECK-LABEL: unsigned_sat_constant_v8i16_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI29_1@toc@ha
; CHECK-NEXT:    vspltisb 5, -1
; CHECK-NEXT:    addi 3, 3, .LCPI29_1@toc@l
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    addis 3, 2, .LCPI29_0@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI29_0@toc@l
; CHECK-NEXT:    vcmpgtuh 3, 2, 3
; CHECK-NEXT:    lvx 4, 0, 3
; CHECK-NEXT:    vadduhm 2, 2, 4
; CHECK-NEXT:    xxsel 34, 34, 37, 35
; CHECK-NEXT:    blr
  %a = add <8 x i16> %x, <i16 42, i16 42, i16 42, i16 42, i16 42, i16 42, i16 42, i16 42>
  %c = icmp ugt <8 x i16> %x, <i16 -43, i16 -43, i16 -43, i16 -43, i16 -43, i16 -43, i16 -43, i16 -43>
  %r = select <8 x i1> %c, <8 x i16> <i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1>, <8 x i16> %a
  ret <8 x i16> %r
}

define <4 x i32> @unsigned_sat_constant_v4i32_using_min(<4 x i32> %x) {
; CHECK-LABEL: unsigned_sat_constant_v4i32_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI30_0@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI30_0@toc@l
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    addis 3, 2, .LCPI30_1@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI30_1@toc@l
; CHECK-NEXT:    vcmpgtuw 4, 3, 2
; CHECK-NEXT:    xxsel 34, 35, 34, 36
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    vadduwm 2, 2, 3
; CHECK-NEXT:    blr
  %c = icmp ult <4 x i32> %x, <i32 -43, i32 -43, i32 -43, i32 -43>
  %s = select <4 x i1> %c, <4 x i32> %x, <4 x i32> <i32 -43, i32 -43, i32 -43, i32 -43>
  %r = add <4 x i32> %s, <i32 42, i32 42, i32 42, i32 42>
  ret <4 x i32> %r
}

define <4 x i32> @unsigned_sat_constant_v4i32_using_cmp_sum(<4 x i32> %x) {
; CHECK-LABEL: unsigned_sat_constant_v4i32_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI31_0@toc@ha
; CHECK-NEXT:    vspltisb 4, -1
; CHECK-NEXT:    addi 3, 3, .LCPI31_0@toc@l
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    vadduwm 3, 2, 3
; CHECK-NEXT:    vcmpgtuw 2, 2, 3
; CHECK-NEXT:    xxsel 34, 35, 36, 34
; CHECK-NEXT:    blr
  %a = add <4 x i32> %x, <i32 42, i32 42, i32 42, i32 42>
  %c = icmp ugt <4 x i32> %x, %a
  %r = select <4 x i1> %c, <4 x i32> <i32 -1, i32 -1, i32 -1, i32 -1>, <4 x i32> %a
  ret <4 x i32> %r
}

define <4 x i32> @unsigned_sat_constant_v4i32_using_cmp_notval(<4 x i32> %x) {
; CHECK-LABEL: unsigned_sat_constant_v4i32_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI32_1@toc@ha
; CHECK-NEXT:    vspltisb 5, -1
; CHECK-NEXT:    addi 3, 3, .LCPI32_1@toc@l
; CHECK-NEXT:    lvx 3, 0, 3
; CHECK-NEXT:    addis 3, 2, .LCPI32_0@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI32_0@toc@l
; CHECK-NEXT:    vcmpgtuw 3, 2, 3
; CHECK-NEXT:    lvx 4, 0, 3
; CHECK-NEXT:    vadduwm 2, 2, 4
; CHECK-NEXT:    xxsel 34, 34, 37, 35
; CHECK-NEXT:    blr
  %a = add <4 x i32> %x, <i32 42, i32 42, i32 42, i32 42>
  %c = icmp ugt <4 x i32> %x, <i32 -43, i32 -43, i32 -43, i32 -43>
  %r = select <4 x i1> %c, <4 x i32> <i32 -1, i32 -1, i32 -1, i32 -1>, <4 x i32> %a
  ret <4 x i32> %r
}

define <2 x i64> @unsigned_sat_constant_v2i64_using_min(<2 x i64> %x) {
; CHECK-LABEL: unsigned_sat_constant_v2i64_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI33_0@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI33_0@toc@l
; CHECK-NEXT:    lxvd2x 0, 0, 3
; CHECK-NEXT:    addis 3, 2, .LCPI33_1@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI33_1@toc@l
; CHECK-NEXT:    xxswapd 35, 0
; CHECK-NEXT:    lxvd2x 0, 0, 3
; CHECK-NEXT:    vcmpgtud 4, 3, 2
; CHECK-NEXT:    xxsel 34, 35, 34, 36
; CHECK-NEXT:    xxswapd 35, 0
; CHECK-NEXT:    vaddudm 2, 2, 3
; CHECK-NEXT:    blr
  %c = icmp ult <2 x i64> %x, <i64 -43, i64 -43>
  %s = select <2 x i1> %c, <2 x i64> %x, <2 x i64> <i64 -43, i64 -43>
  %r = add <2 x i64> %s, <i64 42, i64 42>
  ret <2 x i64> %r
}

define <2 x i64> @unsigned_sat_constant_v2i64_using_cmp_sum(<2 x i64> %x) {
; CHECK-LABEL: unsigned_sat_constant_v2i64_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI34_0@toc@ha
; CHECK-NEXT:    vspltisb 4, -1
; CHECK-NEXT:    addi 3, 3, .LCPI34_0@toc@l
; CHECK-NEXT:    lxvd2x 0, 0, 3
; CHECK-NEXT:    xxswapd 35, 0
; CHECK-NEXT:    vaddudm 3, 2, 3
; CHECK-NEXT:    vcmpgtud 2, 2, 3
; CHECK-NEXT:    xxsel 34, 35, 36, 34
; CHECK-NEXT:    blr
  %a = add <2 x i64> %x, <i64 42, i64 42>
  %c = icmp ugt <2 x i64> %x, %a
  %r = select <2 x i1> %c, <2 x i64> <i64 -1, i64 -1>, <2 x i64> %a
  ret <2 x i64> %r
}

define <2 x i64> @unsigned_sat_constant_v2i64_using_cmp_notval(<2 x i64> %x) {
; CHECK-LABEL: unsigned_sat_constant_v2i64_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    addis 3, 2, .LCPI35_1@toc@ha
; CHECK-NEXT:    vspltisb 5, -1
; CHECK-NEXT:    addi 3, 3, .LCPI35_1@toc@l
; CHECK-NEXT:    lxvd2x 0, 0, 3
; CHECK-NEXT:    addis 3, 2, .LCPI35_0@toc@ha
; CHECK-NEXT:    addi 3, 3, .LCPI35_0@toc@l
; CHECK-NEXT:    xxswapd 35, 0
; CHECK-NEXT:    lxvd2x 0, 0, 3
; CHECK-NEXT:    vcmpgtud 3, 2, 3
; CHECK-NEXT:    xxswapd 36, 0
; CHECK-NEXT:    vaddudm 2, 2, 4
; CHECK-NEXT:    xxsel 34, 34, 37, 35
; CHECK-NEXT:    blr
  %a = add <2 x i64> %x, <i64 42, i64 42>
  %c = icmp ugt <2 x i64> %x, <i64 -43, i64 -43>
  %r = select <2 x i1> %c, <2 x i64> <i64 -1, i64 -1>, <2 x i64> %a
  ret <2 x i64> %r
}

define <16 x i8> @unsigned_sat_variable_v16i8_using_min(<16 x i8> %x, <16 x i8> %y) {
; CHECK-LABEL: unsigned_sat_variable_v16i8_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    xxlnor 36, 35, 35
; CHECK-NEXT:    vcmpgtub 5, 4, 2
; CHECK-NEXT:    xxsel 34, 36, 34, 37
; CHECK-NEXT:    vaddubm 2, 2, 3
; CHECK-NEXT:    blr
  %noty = xor <16 x i8> %y, <i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1>
  %c = icmp ult <16 x i8> %x, %noty
  %s = select <16 x i1> %c, <16 x i8> %x, <16 x i8> %noty
  %r = add <16 x i8> %s, %y
  ret <16 x i8> %r
}

define <16 x i8> @unsigned_sat_variable_v16i8_using_cmp_sum(<16 x i8> %x, <16 x i8> %y) {
; CHECK-LABEL: unsigned_sat_variable_v16i8_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    vaddubm 3, 2, 3
; CHECK-NEXT:    vspltisb 4, -1
; CHECK-NEXT:    vcmpgtub 2, 2, 3
; CHECK-NEXT:    xxsel 34, 35, 36, 34
; CHECK-NEXT:    blr
  %a = add <16 x i8> %x, %y
  %c = icmp ugt <16 x i8> %x, %a
  %r = select <16 x i1> %c, <16 x i8> <i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1>, <16 x i8> %a
  ret <16 x i8> %r
}

define <16 x i8> @unsigned_sat_variable_v16i8_using_cmp_notval(<16 x i8> %x, <16 x i8> %y) {
; CHECK-LABEL: unsigned_sat_variable_v16i8_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    xxlnor 36, 35, 35
; CHECK-NEXT:    vspltisb 5, -1
; CHECK-NEXT:    vcmpgtub 4, 2, 4
; CHECK-NEXT:    vaddubm 2, 2, 3
; CHECK-NEXT:    xxsel 34, 34, 37, 36
; CHECK-NEXT:    blr
  %noty = xor <16 x i8> %y, <i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1>
  %a = add <16 x i8> %x, %y
  %c = icmp ugt <16 x i8> %x, %noty
  %r = select <16 x i1> %c, <16 x i8> <i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1, i8 -1>, <16 x i8> %a
  ret <16 x i8> %r
}

define <8 x i16> @unsigned_sat_variable_v8i16_using_min(<8 x i16> %x, <8 x i16> %y) {
; CHECK-LABEL: unsigned_sat_variable_v8i16_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    xxlnor 36, 35, 35
; CHECK-NEXT:    vcmpgtuh 5, 4, 2
; CHECK-NEXT:    xxsel 34, 36, 34, 37
; CHECK-NEXT:    vadduhm 2, 2, 3
; CHECK-NEXT:    blr
  %noty = xor <8 x i16> %y, <i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1>
  %c = icmp ult <8 x i16> %x, %noty
  %s = select <8 x i1> %c, <8 x i16> %x, <8 x i16> %noty
  %r = add <8 x i16> %s, %y
  ret <8 x i16> %r
}

define <8 x i16> @unsigned_sat_variable_v8i16_using_cmp_sum(<8 x i16> %x, <8 x i16> %y) {
; CHECK-LABEL: unsigned_sat_variable_v8i16_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    vadduhm 3, 2, 3
; CHECK-NEXT:    vspltisb 4, -1
; CHECK-NEXT:    vcmpgtuh 2, 2, 3
; CHECK-NEXT:    xxsel 34, 35, 36, 34
; CHECK-NEXT:    blr
  %a = add <8 x i16> %x, %y
  %c = icmp ugt <8 x i16> %x, %a
  %r = select <8 x i1> %c, <8 x i16> <i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1>, <8 x i16> %a
  ret <8 x i16> %r
}

define <8 x i16> @unsigned_sat_variable_v8i16_using_cmp_notval(<8 x i16> %x, <8 x i16> %y) {
; CHECK-LABEL: unsigned_sat_variable_v8i16_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    xxlnor 36, 35, 35
; CHECK-NEXT:    vspltisb 5, -1
; CHECK-NEXT:    vcmpgtuh 4, 2, 4
; CHECK-NEXT:    vadduhm 2, 2, 3
; CHECK-NEXT:    xxsel 34, 34, 37, 36
; CHECK-NEXT:    blr
  %noty = xor <8 x i16> %y, <i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1>
  %a = add <8 x i16> %x, %y
  %c = icmp ugt <8 x i16> %x, %noty
  %r = select <8 x i1> %c, <8 x i16> <i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1, i16 -1>, <8 x i16> %a
  ret <8 x i16> %r
}

define <4 x i32> @unsigned_sat_variable_v4i32_using_min(<4 x i32> %x, <4 x i32> %y) {
; CHECK-LABEL: unsigned_sat_variable_v4i32_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    xxlnor 36, 35, 35
; CHECK-NEXT:    vcmpgtuw 5, 4, 2
; CHECK-NEXT:    xxsel 34, 36, 34, 37
; CHECK-NEXT:    vadduwm 2, 2, 3
; CHECK-NEXT:    blr
  %noty = xor <4 x i32> %y, <i32 -1, i32 -1, i32 -1, i32 -1>
  %c = icmp ult <4 x i32> %x, %noty
  %s = select <4 x i1> %c, <4 x i32> %x, <4 x i32> %noty
  %r = add <4 x i32> %s, %y
  ret <4 x i32> %r
}

define <4 x i32> @unsigned_sat_variable_v4i32_using_cmp_sum(<4 x i32> %x, <4 x i32> %y) {
; CHECK-LABEL: unsigned_sat_variable_v4i32_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    vadduwm 3, 2, 3
; CHECK-NEXT:    vspltisb 4, -1
; CHECK-NEXT:    vcmpgtuw 2, 2, 3
; CHECK-NEXT:    xxsel 34, 35, 36, 34
; CHECK-NEXT:    blr
  %a = add <4 x i32> %x, %y
  %c = icmp ugt <4 x i32> %x, %a
  %r = select <4 x i1> %c, <4 x i32> <i32 -1, i32 -1, i32 -1, i32 -1>, <4 x i32> %a
  ret <4 x i32> %r
}

define <4 x i32> @unsigned_sat_variable_v4i32_using_cmp_notval(<4 x i32> %x, <4 x i32> %y) {
; CHECK-LABEL: unsigned_sat_variable_v4i32_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    xxlnor 36, 35, 35
; CHECK-NEXT:    vspltisb 5, -1
; CHECK-NEXT:    vcmpgtuw 4, 2, 4
; CHECK-NEXT:    vadduwm 2, 2, 3
; CHECK-NEXT:    xxsel 34, 34, 37, 36
; CHECK-NEXT:    blr
  %noty = xor <4 x i32> %y, <i32 -1, i32 -1, i32 -1, i32 -1>
  %a = add <4 x i32> %x, %y
  %c = icmp ugt <4 x i32> %x, %noty
  %r = select <4 x i1> %c, <4 x i32> <i32 -1, i32 -1, i32 -1, i32 -1>, <4 x i32> %a
  ret <4 x i32> %r
}

define <2 x i64> @unsigned_sat_variable_v2i64_using_min(<2 x i64> %x, <2 x i64> %y) {
; CHECK-LABEL: unsigned_sat_variable_v2i64_using_min:
; CHECK:       # %bb.0:
; CHECK-NEXT:    xxlnor 36, 35, 35
; CHECK-NEXT:    vcmpgtud 5, 4, 2
; CHECK-NEXT:    xxsel 34, 36, 34, 37
; CHECK-NEXT:    vaddudm 2, 2, 3
; CHECK-NEXT:    blr
  %noty = xor <2 x i64> %y, <i64 -1, i64 -1>
  %c = icmp ult <2 x i64> %x, %noty
  %s = select <2 x i1> %c, <2 x i64> %x, <2 x i64> %noty
  %r = add <2 x i64> %s, %y
  ret <2 x i64> %r
}

define <2 x i64> @unsigned_sat_variable_v2i64_using_cmp_sum(<2 x i64> %x, <2 x i64> %y) {
; CHECK-LABEL: unsigned_sat_variable_v2i64_using_cmp_sum:
; CHECK:       # %bb.0:
; CHECK-NEXT:    vaddudm 3, 2, 3
; CHECK-NEXT:    vspltisb 4, -1
; CHECK-NEXT:    vcmpgtud 2, 2, 3
; CHECK-NEXT:    xxsel 34, 35, 36, 34
; CHECK-NEXT:    blr
  %a = add <2 x i64> %x, %y
  %c = icmp ugt <2 x i64> %x, %a
  %r = select <2 x i1> %c, <2 x i64> <i64 -1, i64 -1>, <2 x i64> %a
  ret <2 x i64> %r
}

define <2 x i64> @unsigned_sat_variable_v2i64_using_cmp_notval(<2 x i64> %x, <2 x i64> %y) {
; CHECK-LABEL: unsigned_sat_variable_v2i64_using_cmp_notval:
; CHECK:       # %bb.0:
; CHECK-NEXT:    xxlnor 36, 35, 35
; CHECK-NEXT:    vspltisb 5, -1
; CHECK-NEXT:    vcmpgtud 4, 2, 4
; CHECK-NEXT:    vaddudm 2, 2, 3
; CHECK-NEXT:    xxsel 34, 34, 37, 36
; CHECK-NEXT:    blr
  %noty = xor <2 x i64> %y, <i64 -1, i64 -1>
  %a = add <2 x i64> %x, %y
  %c = icmp ugt <2 x i64> %x, %noty
  %r = select <2 x i1> %c, <2 x i64> <i64 -1, i64 -1>, <2 x i64> %a
  ret <2 x i64> %r
}

