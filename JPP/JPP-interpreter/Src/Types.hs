module Src.Types where

import Src.Jabba.Abs

data VarType
    = VTInt 
    | VTBool 
    | VTString 
    | VTVoid 
    | Fn [FnArg] VarType
    | VTTab VarType
    deriving (Eq)

type FnArg = (VarType, VarMutability, VarRef)


data VarMutability
    = VMMut 
    | VMConst 
    deriving (Eq)


data VarRef
    = VRRef
    | VRCopy
    deriving (Eq)


absTypeToVarType :: Type -> VarType
absTypeToVarType (TInt _) = VTInt
absTypeToVarType (TBool _) = VTBool
absTypeToVarType (TString _) = VTString
absTypeToVarType (TVoid _) = VTVoid
absTypeToVarType (TFun _ args ret) = Fn (map absTypeToFnArg args) (absTypeToVarType ret)
absTypeToVarType (TTab _ t) = VTTab $ absTypeToVarType t


absTypeToFnArg :: TArg -> FnArg
absTypeToFnArg (TRefMutArg    _ t) = (absTypeToVarType t, VMMut,   VRRef)
absTypeToFnArg (TRefConstArg  _ t) = (absTypeToVarType t, VMConst, VRRef)
absTypeToFnArg (TCopyMutArg   _ t) = (absTypeToVarType t, VMMut,   VRCopy)
absTypeToFnArg (TCopyConstArg _ t) = (absTypeToVarType t, VMConst, VRCopy)


instance Show VarType where
    show VTInt = "Integer"
    show VTString = "String"
    show VTBool = "Boolean"
    show VTVoid = "Unit"
    show (Fn args ret) = "(" ++ show args ++ ") -> " ++ show ret
    show (VTTab t) = "[" ++ show t ++ "]"


instance Show VarMutability where
    show VMMut = "mut"
    show VMConst = "const"


instance Show VarRef where
    show VRRef = "ref"
    show VRCopy = "copy"
