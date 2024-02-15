module Src.Errors where

import Src.Jabba.Abs
import Src.Types

type Err = Either ErrHolder

data ErrHolder 
    = ParserErr String
    | TypeChecker BNFC'Position ErrType
    | RuntimeError BNFC'Position RuntimeType
    | ControlledExit Int

data ErrType
    = ImmutVar Ident
    | NotDeclVar Ident
    | NotDeclFun Ident
    | WrongType Ident VarType [VarType]
    | WrongTypeArg Int VarType VarType
    | VarAlreadyDecl Ident
    | WrongNumberOfArgs Int Int
    | ExprRefPass Int
    | ImmutMutPass Int
    | WrongTypeOp String VarType
    | WrongTypeBiOp String VarType VarType
    | TernaryMismatch VarType VarType
    | ZeroLiternalDiv
    | ConstantAssign Ident
    | MismatchedReturnTypes VarType VarType
    | TopLevelProgramReturn VarType
    | TopLevelProgramMaybeReturn VarType
    | ForRangeTypeMismatch VarType VarType
    | TopLevelProgramLoopFlow
    | FunctionLoopFlow Ident
    | FunctionNoReturn Ident VarType
    | FunctionReturnMismatch Ident VarType VarType
    | FunctionMaybeReturn Ident VarType
    | NotFnType VarType
    | FunctionWithoutInitializer Ident
    | InvalidArrayIndex VarType
    | NotAnArray Ident
    | EmptyArray
    | ArrayElementsMismatch [VarType]
    | TableWithoutInitializer Ident
    deriving (Eq)


data RuntimeType
    = ZeroDiv
    | TypeCheckerDidntCatch
    | CannotCast String String
    | AssertionFailed String
    | ArrayIndexOutOfBounds Int Int
    | ArrayInvalidSize Int
    deriving (Eq)


showPos :: BNFC'Position -> String
showPos BNFC'NoPosition = "unknown position"
showPos (BNFC'Position l c) = "line " ++ show l ++ ", column " ++ show c
showPos pos = show pos


showI :: Ident -> String
showI (Ident s) = "\"" ++ s ++ "\""


instance Show ErrHolder where
    show (TypeChecker pos errT) = 
        "There is a problem at " ++ showPos pos ++
        ": " ++ show errT
    show (ParserErr err) = 
        "There was a parsing problem " ++ 
        show err
    show (RuntimeError pos errT) = 
        "There was a runtime problem at " ++ showPos pos ++
        ": " ++ show errT
    show (ControlledExit i) =
        "Program exited with code " ++ show i


instance Show ErrType where
    show (ImmutVar i) = "variable " ++ showI i ++ " is immutable."
    show (NotDeclVar i) = "variable " ++ showI i ++ " was not declared in current scope."
    show (NotDeclFun i) = "function " ++ showI i ++ " was not declared in current scope."
    show (WrongType i t ex) = "variable " ++ showI i ++ " is of type " ++ show t ++ " but expected " ++ show ex ++ "."
    show (WrongTypeArg i t ex) = "argument " ++ show (i + 1) ++ " is of type " ++ show t ++ " but expected " ++ show ex ++ "."
    show (VarAlreadyDecl i) = "variable " ++ showI i ++ " was already declared in declaration block."
    show (WrongNumberOfArgs ex got) = "function was called with " ++ show got ++ " arguments but expected " ++ show ex ++ " arguments."
    show (ExprRefPass i) = "expression " ++ show (i + 1) ++ " was passed to a mutable reference argument."
    show (ImmutMutPass i) = "immutable variable " ++ show (i + 1) ++ " was passed to a mutable reference argument."
    show (WrongTypeOp op t) = "operator " ++ op ++ " is not defined for type " ++ show t ++ "."
    show (WrongTypeBiOp op t1 t2) = "operator " ++ op ++ " is not defined for types " ++ show t1 ++ " and " ++ show t2 ++ "."
    show (TernaryMismatch t1 t2) = "ternary operator has type " ++ show t1 ++ " and " ++ show t2 ++ " as branches."
    show ZeroLiternalDiv = "division by zero literal."
    show (ConstantAssign i) = "variable " ++ showI i ++ " is immutable and cannot be assigned to."
    show (MismatchedReturnTypes t1 t2) = "multiple return statements with types " ++ show t1 ++ " and " ++ show t2 ++ "."
    show (TopLevelProgramReturn t) = "top level program returns type " ++ show t ++ " but expected nothing."
    show (TopLevelProgramMaybeReturn t) = "top level program may return type " ++ show t ++ " but expected nothing."
    show (ForRangeTypeMismatch t1 t2) = "for range has type " ++ show t1 ++ " and " ++ show t2 ++ " as bounds."
    show TopLevelProgramLoopFlow = "top level program contains loop flow statement."
    show (FunctionLoopFlow i) = "function " ++ show i ++ " contains loop flow statement that is not in a loop."
    show (FunctionNoReturn i t) = "function " ++ show i ++ " doens't have a return statement and returns type " ++ show t ++ "."
    show (FunctionReturnMismatch i t1 t2) = "function " ++ show i ++ " returns type " ++ show t1 ++ " but evaluated " ++ show t2 ++ "."
    show (FunctionMaybeReturn i t) = "not all paths in function " ++ show i ++ " return a value and return type is " ++ show t ++ "."
    show (NotFnType t) = "type " ++ show t ++ " is not a function type."
    show (FunctionWithoutInitializer i) = "function " ++ show i ++ " is not initialized."
    show (InvalidArrayIndex t) = "array index is of type " ++ show t ++ " but expected Integer."
    show (NotAnArray i) = "variable " ++ show i ++ " is not an array."
    show EmptyArray = "array is empty - cannot deduce type."
    show (ArrayElementsMismatch ts) = "array elements have types " ++ show ts ++ " but expected all to be the same."
    show (TableWithoutInitializer i) = "table " ++ show i ++ " is not initialized."


instance Show RuntimeType where
    show ZeroDiv = "division by zero."
    show TypeCheckerDidntCatch = "type checker didn't catch this error."
    show (CannotCast s t) = "cannot cast " ++ s ++ " to type " ++ show t ++ "."
    show (AssertionFailed s) = "assertion failed:\n" ++ s ++ "."
    show (ArrayIndexOutOfBounds i l) = "array index " ++ show i ++ " is out of bounds (length " ++ show l ++ ")."
    show (ArrayInvalidSize i) = "array size " ++ show i ++ " is invalid."
    