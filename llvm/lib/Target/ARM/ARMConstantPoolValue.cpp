//===- ARMConstantPoolValue.cpp - ARM constantpool value --------*- C++ -*-===//
//
//                     The LLVM Compiler Infrastructure
//
// This file was developed by Evan Cheng and is distributed under the
// University of Illinois Open Source License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// This file implements the ARM specific constantpool value class.
//
//===----------------------------------------------------------------------===//

#include "ARMConstantPoolValue.h"
#include "llvm/ADT/FoldingSet.h"
#include "llvm/GlobalValue.h"
#include "llvm/Type.h"
using namespace llvm;

ARMConstantPoolValue::ARMConstantPoolValue(GlobalValue *gv, unsigned id,
                                           ARMCP::ARMCPKind k,
                                           unsigned char PCAdj)
  : MachineConstantPoolValue((const Type*)gv->getType()),
    GV(gv), S(NULL), LabelId(id), Kind(k), PCAdjust(PCAdj) {}

ARMConstantPoolValue::ARMConstantPoolValue(const char *s, unsigned id,
                                           ARMCP::ARMCPKind k,
                                           unsigned char PCAdj)
  : MachineConstantPoolValue((const Type*)Type::Int32Ty),
    GV(NULL), S(s), LabelId(id), Kind(k), PCAdjust(PCAdj) {}

int ARMConstantPoolValue::getExistingMachineCPValue(MachineConstantPool *CP,
                                                    unsigned Alignment) {
  unsigned AlignMask = (1 << Alignment)-1;
  const std::vector<MachineConstantPoolEntry> Constants = CP->getConstants();
  for (unsigned i = 0, e = Constants.size(); i != e; ++i) {
    if (Constants[i].isMachineConstantPoolEntry() &&
        (Constants[i].Offset & AlignMask) == 0) {
      ARMConstantPoolValue *CPV =
        (ARMConstantPoolValue *)Constants[i].Val.MachineCPVal;
      if (CPV->GV == GV &&
          CPV->S == S &&
          CPV->LabelId == LabelId &&
          CPV->Kind == Kind &&
          CPV->PCAdjust == PCAdjust)
        return i;
    }
  }

  return -1;
}

void
ARMConstantPoolValue::AddSelectionDAGCSEId(FoldingSetNodeID &ID) {
  ID.AddPointer(GV);
  ID.AddPointer(S);
  ID.AddInteger(LabelId);
  ID.AddInteger((unsigned)Kind);
  ID.AddInteger(PCAdjust);
}

void ARMConstantPoolValue::print(std::ostream &O) const {
  if (GV)
    O << GV->getName();
  else
    O << S;
  if (isNonLazyPointer()) O << "$non_lazy_ptr";
  else if (isStub()) O << "$stub";
  if (PCAdjust != 0) O << "-(LPIC" << LabelId << "+"
                       << (unsigned)PCAdjust << ")";
}
