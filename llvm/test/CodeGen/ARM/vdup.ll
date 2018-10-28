; NOTE: Assertions have been autogenerated by utils/update_llc_test_checks.py
; RUN: llc < %s -mtriple=arm-eabi -float-abi=soft -mattr=+neon -verify-machineinstrs | FileCheck %s

define <8 x i8> @v_dup8(i8 %A) nounwind {
; CHECK-LABEL: v_dup8:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.8 d16, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <8 x i8> zeroinitializer, i8 %A, i32 0
	%tmp2 = insertelement <8 x i8> %tmp1, i8 %A, i32 1
	%tmp3 = insertelement <8 x i8> %tmp2, i8 %A, i32 2
	%tmp4 = insertelement <8 x i8> %tmp3, i8 %A, i32 3
	%tmp5 = insertelement <8 x i8> %tmp4, i8 %A, i32 4
	%tmp6 = insertelement <8 x i8> %tmp5, i8 %A, i32 5
	%tmp7 = insertelement <8 x i8> %tmp6, i8 %A, i32 6
	%tmp8 = insertelement <8 x i8> %tmp7, i8 %A, i32 7
	ret <8 x i8> %tmp8
}

define <4 x i16> @v_dup16(i16 %A) nounwind {
; CHECK-LABEL: v_dup16:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.16 d16, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <4 x i16> zeroinitializer, i16 %A, i32 0
	%tmp2 = insertelement <4 x i16> %tmp1, i16 %A, i32 1
	%tmp3 = insertelement <4 x i16> %tmp2, i16 %A, i32 2
	%tmp4 = insertelement <4 x i16> %tmp3, i16 %A, i32 3
	ret <4 x i16> %tmp4
}

define <2 x i32> @v_dup32(i32 %A) nounwind {
; CHECK-LABEL: v_dup32:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.32 d16, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <2 x i32> zeroinitializer, i32 %A, i32 0
	%tmp2 = insertelement <2 x i32> %tmp1, i32 %A, i32 1
	ret <2 x i32> %tmp2
}

define <2 x float> @v_dupfloat(float %A) nounwind {
; CHECK-LABEL: v_dupfloat:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.32 d16, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <2 x float> zeroinitializer, float %A, i32 0
	%tmp2 = insertelement <2 x float> %tmp1, float %A, i32 1
	ret <2 x float> %tmp2
}

define <16 x i8> @v_dupQ8(i8 %A) nounwind {
; CHECK-LABEL: v_dupQ8:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.8 q8, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <16 x i8> zeroinitializer, i8 %A, i32 0
	%tmp2 = insertelement <16 x i8> %tmp1, i8 %A, i32 1
	%tmp3 = insertelement <16 x i8> %tmp2, i8 %A, i32 2
	%tmp4 = insertelement <16 x i8> %tmp3, i8 %A, i32 3
	%tmp5 = insertelement <16 x i8> %tmp4, i8 %A, i32 4
	%tmp6 = insertelement <16 x i8> %tmp5, i8 %A, i32 5
	%tmp7 = insertelement <16 x i8> %tmp6, i8 %A, i32 6
	%tmp8 = insertelement <16 x i8> %tmp7, i8 %A, i32 7
	%tmp9 = insertelement <16 x i8> %tmp8, i8 %A, i32 8
	%tmp10 = insertelement <16 x i8> %tmp9, i8 %A, i32 9
	%tmp11 = insertelement <16 x i8> %tmp10, i8 %A, i32 10
	%tmp12 = insertelement <16 x i8> %tmp11, i8 %A, i32 11
	%tmp13 = insertelement <16 x i8> %tmp12, i8 %A, i32 12
	%tmp14 = insertelement <16 x i8> %tmp13, i8 %A, i32 13
	%tmp15 = insertelement <16 x i8> %tmp14, i8 %A, i32 14
	%tmp16 = insertelement <16 x i8> %tmp15, i8 %A, i32 15
	ret <16 x i8> %tmp16
}

define <8 x i16> @v_dupQ16(i16 %A) nounwind {
; CHECK-LABEL: v_dupQ16:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.16 q8, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <8 x i16> zeroinitializer, i16 %A, i32 0
	%tmp2 = insertelement <8 x i16> %tmp1, i16 %A, i32 1
	%tmp3 = insertelement <8 x i16> %tmp2, i16 %A, i32 2
	%tmp4 = insertelement <8 x i16> %tmp3, i16 %A, i32 3
	%tmp5 = insertelement <8 x i16> %tmp4, i16 %A, i32 4
	%tmp6 = insertelement <8 x i16> %tmp5, i16 %A, i32 5
	%tmp7 = insertelement <8 x i16> %tmp6, i16 %A, i32 6
	%tmp8 = insertelement <8 x i16> %tmp7, i16 %A, i32 7
	ret <8 x i16> %tmp8
}

define <4 x i32> @v_dupQ32(i32 %A) nounwind {
; CHECK-LABEL: v_dupQ32:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.32 q8, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <4 x i32> zeroinitializer, i32 %A, i32 0
	%tmp2 = insertelement <4 x i32> %tmp1, i32 %A, i32 1
	%tmp3 = insertelement <4 x i32> %tmp2, i32 %A, i32 2
	%tmp4 = insertelement <4 x i32> %tmp3, i32 %A, i32 3
	ret <4 x i32> %tmp4
}

define <4 x float> @v_dupQfloat(float %A) nounwind {
; CHECK-LABEL: v_dupQfloat:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.32 q8, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <4 x float> zeroinitializer, float %A, i32 0
	%tmp2 = insertelement <4 x float> %tmp1, float %A, i32 1
	%tmp3 = insertelement <4 x float> %tmp2, float %A, i32 2
	%tmp4 = insertelement <4 x float> %tmp3, float %A, i32 3
	ret <4 x float> %tmp4
}

; Check to make sure it works with shuffles, too.

define <8 x i8> @v_shuffledup8(i8 %A) nounwind {
; CHECK-LABEL: v_shuffledup8:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.8 d16, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <8 x i8> undef, i8 %A, i32 0
	%tmp2 = shufflevector <8 x i8> %tmp1, <8 x i8> undef, <8 x i32> zeroinitializer
	ret <8 x i8> %tmp2
}

define <4 x i16> @v_shuffledup16(i16 %A) nounwind {
; CHECK-LABEL: v_shuffledup16:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.16 d16, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <4 x i16> undef, i16 %A, i32 0
	%tmp2 = shufflevector <4 x i16> %tmp1, <4 x i16> undef, <4 x i32> zeroinitializer
	ret <4 x i16> %tmp2
}

define <2 x i32> @v_shuffledup32(i32 %A) nounwind {
; CHECK-LABEL: v_shuffledup32:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.32 d16, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <2 x i32> undef, i32 %A, i32 0
	%tmp2 = shufflevector <2 x i32> %tmp1, <2 x i32> undef, <2 x i32> zeroinitializer
	ret <2 x i32> %tmp2
}

define <2 x float> @v_shuffledupfloat(float %A) nounwind {
; CHECK-LABEL: v_shuffledupfloat:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.32 d16, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <2 x float> undef, float %A, i32 0
	%tmp2 = shufflevector <2 x float> %tmp1, <2 x float> undef, <2 x i32> zeroinitializer
	ret <2 x float> %tmp2
}

define <16 x i8> @v_shuffledupQ8(i8 %A) nounwind {
; CHECK-LABEL: v_shuffledupQ8:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.8 q8, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <16 x i8> undef, i8 %A, i32 0
	%tmp2 = shufflevector <16 x i8> %tmp1, <16 x i8> undef, <16 x i32> zeroinitializer
	ret <16 x i8> %tmp2
}

define <8 x i16> @v_shuffledupQ16(i16 %A) nounwind {
; CHECK-LABEL: v_shuffledupQ16:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.16 q8, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <8 x i16> undef, i16 %A, i32 0
	%tmp2 = shufflevector <8 x i16> %tmp1, <8 x i16> undef, <8 x i32> zeroinitializer
	ret <8 x i16> %tmp2
}

define <4 x i32> @v_shuffledupQ32(i32 %A) nounwind {
; CHECK-LABEL: v_shuffledupQ32:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.32 q8, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <4 x i32> undef, i32 %A, i32 0
	%tmp2 = shufflevector <4 x i32> %tmp1, <4 x i32> undef, <4 x i32> zeroinitializer
	ret <4 x i32> %tmp2
}

define <4 x float> @v_shuffledupQfloat(float %A) nounwind {
; CHECK-LABEL: v_shuffledupQfloat:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.32 q8, r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = insertelement <4 x float> undef, float %A, i32 0
	%tmp2 = shufflevector <4 x float> %tmp1, <4 x float> undef, <4 x i32> zeroinitializer
	ret <4 x float> %tmp2
}

define <8 x i8> @vduplane8(<8 x i8>* %A) nounwind {
; CHECK-LABEL: vduplane8:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vldr d16, [r0]
; CHECK-NEXT:    vdup.8 d16, d16[1]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = load <8 x i8>, <8 x i8>* %A
	%tmp2 = shufflevector <8 x i8> %tmp1, <8 x i8> undef, <8 x i32> < i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1 >
	ret <8 x i8> %tmp2
}

define <4 x i16> @vduplane16(<4 x i16>* %A) nounwind {
; CHECK-LABEL: vduplane16:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vldr d16, [r0]
; CHECK-NEXT:    vdup.16 d16, d16[1]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = load <4 x i16>, <4 x i16>* %A
	%tmp2 = shufflevector <4 x i16> %tmp1, <4 x i16> undef, <4 x i32> < i32 1, i32 1, i32 1, i32 1 >
	ret <4 x i16> %tmp2
}

define <2 x i32> @vduplane32(<2 x i32>* %A) nounwind {
; CHECK-LABEL: vduplane32:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vldr d16, [r0]
; CHECK-NEXT:    vdup.32 d16, d16[1]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = load <2 x i32>, <2 x i32>* %A
	%tmp2 = shufflevector <2 x i32> %tmp1, <2 x i32> undef, <2 x i32> < i32 1, i32 1 >
	ret <2 x i32> %tmp2
}

define <2 x float> @vduplanefloat(<2 x float>* %A) nounwind {
; CHECK-LABEL: vduplanefloat:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vldr d16, [r0]
; CHECK-NEXT:    vdup.32 d16, d16[1]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
	%tmp1 = load <2 x float>, <2 x float>* %A
	%tmp2 = shufflevector <2 x float> %tmp1, <2 x float> undef, <2 x i32> < i32 1, i32 1 >
	ret <2 x float> %tmp2
}

define <16 x i8> @vduplaneQ8(<8 x i8>* %A) nounwind {
; CHECK-LABEL: vduplaneQ8:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vldr d16, [r0]
; CHECK-NEXT:    vdup.8 q8, d16[1]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = load <8 x i8>, <8 x i8>* %A
	%tmp2 = shufflevector <8 x i8> %tmp1, <8 x i8> undef, <16 x i32> < i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1 >
	ret <16 x i8> %tmp2
}

define <8 x i16> @vduplaneQ16(<4 x i16>* %A) nounwind {
; CHECK-LABEL: vduplaneQ16:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vldr d16, [r0]
; CHECK-NEXT:    vdup.16 q8, d16[1]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = load <4 x i16>, <4 x i16>* %A
	%tmp2 = shufflevector <4 x i16> %tmp1, <4 x i16> undef, <8 x i32> < i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1, i32 1 >
	ret <8 x i16> %tmp2
}

define <4 x i32> @vduplaneQ32(<2 x i32>* %A) nounwind {
; CHECK-LABEL: vduplaneQ32:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vldr d16, [r0]
; CHECK-NEXT:    vdup.32 q8, d16[1]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = load <2 x i32>, <2 x i32>* %A
	%tmp2 = shufflevector <2 x i32> %tmp1, <2 x i32> undef, <4 x i32> < i32 1, i32 1, i32 1, i32 1 >
	ret <4 x i32> %tmp2
}

define <4 x float> @vduplaneQfloat(<2 x float>* %A) nounwind {
; CHECK-LABEL: vduplaneQfloat:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vldr d16, [r0]
; CHECK-NEXT:    vdup.32 q8, d16[1]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
	%tmp1 = load <2 x float>, <2 x float>* %A
	%tmp2 = shufflevector <2 x float> %tmp1, <2 x float> undef, <4 x i32> < i32 1, i32 1, i32 1, i32 1 >
	ret <4 x float> %tmp2
}

define <2 x i64> @foo(<2 x i64> %arg0_int64x1_t) nounwind readnone {
; CHECK-LABEL: foo:
; CHECK:       @ %bb.0: @ %entry
; CHECK-NEXT:    mov r0, r2
; CHECK-NEXT:    mov r1, r3
; CHECK-NEXT:    mov pc, lr
entry:
  %0 = shufflevector <2 x i64> %arg0_int64x1_t, <2 x i64> undef, <2 x i32> <i32 1, i32 1>
  ret <2 x i64> %0
}

define <2 x i64> @bar(<2 x i64> %arg0_int64x1_t) nounwind readnone {
; CHECK-LABEL: bar:
; CHECK:       @ %bb.0: @ %entry
; CHECK-NEXT:    mov r2, r0
; CHECK-NEXT:    mov r3, r1
; CHECK-NEXT:    mov pc, lr
entry:
  %0 = shufflevector <2 x i64> %arg0_int64x1_t, <2 x i64> undef, <2 x i32> <i32 0, i32 0>
  ret <2 x i64> %0
}

define <2 x double> @baz(<2 x double> %arg0_int64x1_t) nounwind readnone {
; CHECK-LABEL: baz:
; CHECK:       @ %bb.0: @ %entry
; CHECK-NEXT:    mov r0, r2
; CHECK-NEXT:    mov r1, r3
; CHECK-NEXT:    mov pc, lr
entry:
  %0 = shufflevector <2 x double> %arg0_int64x1_t, <2 x double> undef, <2 x i32> <i32 1, i32 1>
  ret <2 x double> %0
}

define <2 x double> @qux(<2 x double> %arg0_int64x1_t) nounwind readnone {
; CHECK-LABEL: qux:
; CHECK:       @ %bb.0: @ %entry
; CHECK-NEXT:    mov r2, r0
; CHECK-NEXT:    mov r3, r1
; CHECK-NEXT:    mov pc, lr
entry:
  %0 = shufflevector <2 x double> %arg0_int64x1_t, <2 x double> undef, <2 x i32> <i32 0, i32 0>
  ret <2 x double> %0
}

; Radar 7373643
define void @redundantVdup(<8 x i8>* %ptr) nounwind {
; CHECK-LABEL: redundantVdup:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vmov.i8 d16, #0x80
; CHECK-NEXT:    vstr d16, [r0]
; CHECK-NEXT:    mov pc, lr
  %1 = insertelement <8 x i8> undef, i8 -128, i32 0
  %2 = shufflevector <8 x i8> %1, <8 x i8> undef, <8 x i32> zeroinitializer
  store <8 x i8> %2, <8 x i8>* %ptr, align 8
  ret void
}

define <4 x i32> @tdupi(i32 %x, i32 %y) {
; CHECK-LABEL: tdupi:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.32 q8, r0
; CHECK-NEXT:    vmov.32 d17[1], r1
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
  %1 = insertelement <4 x i32> undef, i32 %x, i32 0
  %2 = insertelement <4 x i32> %1, i32 %x, i32 1
  %3 = insertelement <4 x i32> %2, i32 %x, i32 2
  %4 = insertelement <4 x i32> %3, i32 %y, i32 3
  ret <4 x i32> %4
}

define <4 x float> @tdupf(float %x, float %y) {
; CHECK-LABEL: tdupf:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vdup.32 q0, r0
; CHECK-NEXT:    vmov s3, r1
; CHECK-NEXT:    vmov r0, r1, d0
; CHECK-NEXT:    vmov r2, r3, d1
; CHECK-NEXT:    mov pc, lr
  %1 = insertelement <4 x float> undef, float %x, i32 0
  %2 = insertelement <4 x float> %1, float %x, i32 1
  %3 = insertelement <4 x float> %2, float %x, i32 2
  %4 = insertelement <4 x float> %3, float %y, i32 3
  ret <4 x float> %4
}

; This test checks that when splatting an element from a vector into another,
; the value isn't moved out to GPRs first.
define <4 x i32> @tduplane(<4 x i32> %invec) {
; CHECK-LABEL: tduplane:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vmov d16, r0, r1
; CHECK-NEXT:    mov r0, #255
; CHECK-NEXT:    vdup.32 q8, d16[1]
; CHECK-NEXT:    vmov.32 d17[1], r0
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
  %in = extractelement <4 x i32> %invec, i32 1
  %1 = insertelement <4 x i32> undef, i32 %in, i32 0
  %2 = insertelement <4 x i32> %1, i32 %in, i32 1
  %3 = insertelement <4 x i32> %2, i32 %in, i32 2
  %4 = insertelement <4 x i32> %3, i32 255, i32 3
  ret <4 x i32> %4
}

define <2 x float> @check_f32(<4 x float> %v) nounwind {
; CHECK-LABEL: check_f32:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vmov d17, r2, r3
; CHECK-NEXT:    vmov d16, r0, r1
; CHECK-NEXT:    vdup.32 d16, d17[1]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
  %x = extractelement <4 x float> %v, i32 3
  %1 = insertelement  <2 x float> undef, float %x, i32 0
  %2 = insertelement  <2 x float> %1, float %x, i32 1
  ret <2 x float> %2
}

define <2 x i32> @check_i32(<4 x i32> %v) nounwind {
; CHECK-LABEL: check_i32:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vmov d17, r2, r3
; CHECK-NEXT:    vmov d16, r0, r1
; CHECK-NEXT:    vdup.32 d16, d17[1]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
  %x = extractelement <4 x i32> %v, i32 3
  %1 = insertelement  <2 x i32> undef, i32 %x, i32 0
  %2 = insertelement  <2 x i32> %1, i32 %x, i32 1
  ret <2 x i32> %2
}

define <4 x i16> @check_i16(<8 x i16> %v) nounwind {
; CHECK-LABEL: check_i16:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vmov d17, r2, r3
; CHECK-NEXT:    vmov d16, r0, r1
; CHECK-NEXT:    vdup.16 d16, d16[3]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
  %x = extractelement <8 x i16> %v, i32 3
  %1 = insertelement  <4 x i16> undef, i16 %x, i32 0
  %2 = insertelement  <4 x i16> %1, i16 %x, i32 1
  ret <4 x i16> %2
}

define <8 x i8> @check_i8(<16 x i8> %v) nounwind {
; CHECK-LABEL: check_i8:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    vmov d17, r2, r3
; CHECK-NEXT:    vmov d16, r0, r1
; CHECK-NEXT:    vdup.8 d16, d16[3]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
  %x = extractelement <16 x i8> %v, i32 3
  %1 = insertelement  <8  x i8> undef, i8 %x, i32 0
  %2 = insertelement  <8  x i8> %1, i8 %x, i32 1
  ret <8 x i8> %2
}

; Check that an SPR splat produces a vdup.

define <2 x float> @check_spr_splat2(<2 x float> %p, i16 %q) {
; CHECK-LABEL: check_spr_splat2:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    lsl r2, r2, #16
; CHECK-NEXT:    vmov d17, r0, r1
; CHECK-NEXT:    asr r2, r2, #16
; CHECK-NEXT:    vdup.32 d16, r2
; CHECK-NEXT:    vcvt.f32.s32 d16, d16
; CHECK-NEXT:    vsub.f32 d16, d16, d17
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov pc, lr
  %conv = sitofp i16 %q to float
  %splat.splatinsert = insertelement <2 x float> undef, float %conv, i32 0
  %splat.splat = shufflevector <2 x float> %splat.splatinsert, <2 x float> undef, <2 x i32> zeroinitializer
  %sub = fsub <2 x float> %splat.splat, %p
  ret <2 x float> %sub
}

define <4 x float> @check_spr_splat4(<4 x float> %p, i16 %q) {
; CHECK-LABEL: check_spr_splat4:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    mov r12, sp
; CHECK-NEXT:    vmov d19, r2, r3
; CHECK-NEXT:    vld1.16 {d16[]}, [r12:16]
; CHECK-NEXT:    vmov d18, r0, r1
; CHECK-NEXT:    vmovl.s16 q8, d16
; CHECK-NEXT:    vcvt.f32.s32 q8, q8
; CHECK-NEXT:    vsub.f32 q8, q8, q9
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
  %conv = sitofp i16 %q to float
  %splat.splatinsert = insertelement <4 x float> undef, float %conv, i32 0
  %splat.splat = shufflevector <4 x float> %splat.splatinsert, <4 x float> undef, <4 x i32> zeroinitializer
  %sub = fsub <4 x float> %splat.splat, %p
  ret <4 x float> %sub
}
; Same codegen as above test; scalar is splatted using vld1, so shuffle index is irrelevant.
define <4 x float> @check_spr_splat4_lane1(<4 x float> %p, i16 %q) {
; CHECK-LABEL: check_spr_splat4_lane1:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    mov r12, sp
; CHECK-NEXT:    vmov d19, r2, r3
; CHECK-NEXT:    vld1.16 {d16[]}, [r12:16]
; CHECK-NEXT:    vmov d18, r0, r1
; CHECK-NEXT:    vmovl.s16 q8, d16
; CHECK-NEXT:    vcvt.f32.s32 q8, q8
; CHECK-NEXT:    vsub.f32 q8, q8, q9
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    vmov r2, r3, d17
; CHECK-NEXT:    mov pc, lr
  %conv = sitofp i16 %q to float
  %splat.splatinsert = insertelement <4 x float> undef, float %conv, i32 1
  %splat.splat = shufflevector <4 x float> %splat.splatinsert, <4 x float> undef, <4 x i32> <i32 1, i32 1, i32 1, i32 1>
  %sub = fsub <4 x float> %splat.splat, %p
  ret <4 x float> %sub
}

; Also make sure we don't barf on variable-index extractelts, where we almost
; could have generated a vdup.

define <8 x i8> @check_i8_varidx(<16 x i8> %v, i32 %idx) {
; CHECK-LABEL: check_i8_varidx:
; CHECK:       @ %bb.0:
; CHECK-NEXT:    .save {r11}
; CHECK-NEXT:    push {r11}
; CHECK-NEXT:    .setfp r11, sp
; CHECK-NEXT:    mov r11, sp
; CHECK-NEXT:    .pad #28
; CHECK-NEXT:    sub sp, sp, #28
; CHECK-NEXT:    bic sp, sp, #15
; CHECK-NEXT:    ldr r12, [r11, #4]
; CHECK-NEXT:    vmov d17, r2, r3
; CHECK-NEXT:    vmov d16, r0, r1
; CHECK-NEXT:    mov r1, sp
; CHECK-NEXT:    and r0, r12, #15
; CHECK-NEXT:    vst1.64 {d16, d17}, [r1:128], r0
; CHECK-NEXT:    vld1.8 {d16[]}, [r1]
; CHECK-NEXT:    vmov r0, r1, d16
; CHECK-NEXT:    mov sp, r11
; CHECK-NEXT:    pop {r11}
; CHECK-NEXT:    mov pc, lr
  %x = extractelement <16 x i8> %v, i32 %idx
  %1 = insertelement  <8 x i8> undef, i8 %x, i32 0
  %2 = insertelement  <8 x i8> %1, i8 %x, i32 1
  ret <8 x i8> %2
}
