module Frontend.TypeErrors where

import Data.List
import Frontend.Commons

data ErrHolder
  = DuplicateArgName Pos Pos String
  | NoMain
  | InvalidMainArgs Pos
  | InvalidMainRetType Pos
  | CycleInheritance String Pos
  | UnknownClass String Pos
  | UnknownParentClass String Pos
  | DuplicateClassFieldName String Pos Pos
  | DuplicateClassMethodName String Pos Pos
  | DuplicateSuperClassFieldName String String [(String, Pos, Pos)]
  | InvalidVoidType Pos
  | InvalidFunType Pos
  | MissingReturnStatementInFunctionBlock Pos String
  | MissingReturnStatementInFunctionBlockSometimes Pos String
  | ReturnTypeMismatch String Pos Pos
  | ReturnTypeMismatchKnown String Pos String String
  | DuplicateLocalVarName String [Pos]
  | ConditionalMistype Pos Pos
  | UnknownVar String Pos
  | UnknownFunction String Pos
  | InvalidFunctionArgType Pos String String String
  | InvalidFunctionArgCount Pos String
  | InvalidArraySize Pos String
  | InvalidArrayIndex Pos String
  | InvalidArrayAccess Pos String
  | InvalidClassFieldAccess Pos String
  | UnknownClassField String String Pos
  | InvalidClassMethodAccess Pos String
  | UnknownClassMethod String String Pos
  | InvalidNegType Pos String
  | InvalidNotType Pos String
  | InvalidMulType Pos String String
  | InvalidAddType Pos String String
  | InvalidRelType Pos String String
  | InvalidRelOp Pos String String
  | InvalidRelTypeVoid Pos String String
  | InvalidAndType Pos String String
  | InvalidOrType Pos String String
  | CannotAssignType Pos String String
  | CannotIncrementNotInt Pos String
  | CannotDecrementNotInt Pos String
  | IfCondShouldBeBoolean Pos String
  | WhileCondShouldBeBoolean Pos String
  | ForIterableNotArray Pos String
  | InvalidNullCast Pos String
  | CannotChangeNotLValue Pos
  | IntOverflow Pos Integer
  | ForIterableTypeMismatch Pos String String
  | DeclTypeMismatch Pos String String
  | InvalidFieldName Pos String
  | DuplicateFunctionName String [Pos]
  | DuplicateClassName String [Pos]
  | DuplicateBuiltinFunctionName String [Pos]
  | ReturnVoid Pos
  | InvalidReturnType Pos String String
  | DivisionByZero Pos
  | InvalidSuperClassMethodRetType String Pos String Pos String
  | InvalidSuperClassMethodArgs String Pos [String] Pos [String]
  | InvalidClassName Pos String
  | InvalidArrayType Pos String
  deriving (Eq, Ord, Show)

printErrHolder :: String -> ErrHolder -> IO ()
printErrHolder code (DuplicateArgName firstPos nextPos argName) = do
  putStrLnStderr $ "Error: Duplicate argument name: " ++ argName ++ " at " ++ posToString firstPos ++ " and " ++ posToString nextPos
  printCode code [firstPos, nextPos]
printErrHolder code NoMain = do
  putStrLnStderr "Error: No main function"
  printCode code []
printErrHolder code (InvalidMainArgs pos) = do
  putStrLnStderr $ "Error: Invalid main function arguments at " ++ posToString pos
  putStrLnStderr "Function should have no arguments"
  printCode code [pos]
printErrHolder code (InvalidMainRetType pos) = do
  putStrLnStderr $ "Error: Invalid main function return type at " ++ posToString pos
  putStrLnStderr "Function should return int"
  printCode code [pos]
printErrHolder code (CycleInheritance className pos) = do
  putStrLnStderr $ "Error: Cycle in inheritance of class " ++ className ++ " at " ++ posToString pos
  printCode code [pos]
printErrHolder code (UnknownClass className pos) = do
  putStrLnStderr $ "Error: Unknown class " ++ className ++ " at " ++ posToString pos
  printCode code [pos]
printErrHolder code (UnknownParentClass className pos) = do
  putStrLnStderr $ "Error: Unknown parent class " ++ className ++ " at " ++ posToString pos
  printCode code [pos]
printErrHolder code (DuplicateClassFieldName fieldName firstPos nextPos) = do
  putStrLnStderr $ "Error: Duplicate class field name: " ++ fieldName ++ " at " ++ posToString firstPos ++ " and " ++ posToString nextPos
  printCode code [firstPos, nextPos]
printErrHolder code (DuplicateClassMethodName methodName firstPos nextPos) = do
  putStrLnStderr $ "Error: Duplicate class method name: " ++ methodName ++ " at " ++ posToString firstPos ++ " and " ++ posToString nextPos
  printCode code [firstPos, nextPos]
printErrHolder code (DuplicateSuperClassFieldName className superClassName positions) = do
  mapM_ printSingleVar positions
  where
    printSingleVar :: (String, Pos, Pos) -> IO ()
    printSingleVar (name, p1, p2) = do
      putStrLnStderr ""
      putStrLnStderr $ "Error: Duplicate class field name: " ++ name ++ " in class " ++ className ++ " with super class " ++ superClassName
      putStrLnStderr $ "  " ++ name ++ " at " ++ posToString p1 ++ " and " ++ posToString p2
      printCode code [p1, p2]
printErrHolder code (InvalidVoidType pos) = do
  putStrLnStderr $ "Error: Invalid void type at " ++ posToString pos
  putStrLnStderr "Void type can be used only as a return type of a function"
  printCode code [pos]
printErrHolder code (InvalidFunType pos) = do
  putStrLnStderr $ "Error: Invalid function type at " ++ posToString pos
  putStrLnStderr "Function type cannot be used as a type of a variable"
  printCode code [pos]
printErrHolder code (MissingReturnStatementInFunctionBlock pos funName) = do
  putStrLnStderr $ "Error: Missing return statement in function block at " ++ posToString pos
  putStrLnStderr $ "Function " ++ funName ++ " should return a value"
  printCode code [pos]
printErrHolder code (MissingReturnStatementInFunctionBlockSometimes pos funName) = do
  putStrLnStderr $ "Error: Missing return statement in all branches of function block at " ++ posToString pos
  putStrLnStderr $ "Function " ++ funName ++ " should ALWAYS return a value"
  printCode code [pos]
printErrHolder code (ReturnTypeMismatch funName pos1 pos2) = do
  putStrLnStderr $ "Error: Return type mismatch at " ++ posToString pos1 ++ " and " ++ posToString pos2
  putStrLnStderr $ "Function " ++ funName ++ " returns different types"
  printCode code [pos1, pos2]
printErrHolder code (ReturnTypeMismatchKnown funName pos expectedType actualType) = do
  putStrLnStderr $ "Error: Return type mismatch at " ++ posToString pos
  putStrLnStderr $ "Function " ++ funName ++ " should return a value of type " ++ expectedType ++ " but returns " ++ actualType
  printCode code [pos]
printErrHolder code (DuplicateLocalVarName varName poses) = do
  putStrLnStderr $ "Error: Duplicate local variable name: " ++ varName
  mapM_ (putStrLnStderr . ("  " ++) . posToString) poses
  printCode code poses
printErrHolder code (ConditionalMistype pos1 pos2) = do
  putStrLnStderr $ "Error: Conditional mistype at " ++ posToString pos1 ++ " and " ++ posToString pos2
  printCode code [pos1, pos2]
printErrHolder code (UnknownVar varName pos) = do
  putStrLnStderr $ "Error: Unknown variable " ++ varName ++ " at " ++ posToString pos
  printCode code [pos]
printErrHolder code (UnknownFunction funName pos) = do
  putStrLnStderr $ "Error: Unknown function " ++ funName ++ " at " ++ posToString pos
  printCode code [pos]
printErrHolder code (InvalidFunctionArgType pos funName expectedType actualType) = do
  putStrLnStderr $ "Error: Invalid argument type at " ++ posToString pos
  putStrLnStderr $ "Function " ++ funName ++ " expects argument of type " ++ expectedType ++ " but got " ++ actualType
  printCode code [pos]
printErrHolder code (InvalidFunctionArgCount pos funName) = do
  putStrLnStderr $ "Error: Invalid argument count at " ++ posToString pos
  putStrLnStderr $ "Function " ++ funName ++ " expects different number of arguments"
  printCode code [pos]
printErrHolder code (InvalidArraySize pos actualType) = do
  putStrLnStderr $ "Error: Invalid array size at " ++ posToString pos
  putStrLnStderr $ "Array size should be of type int but got " ++ actualType
  printCode code [pos]
printErrHolder code (InvalidArrayIndex pos actualType) = do
  putStrLnStderr $ "Error: Invalid array index at " ++ posToString pos
  putStrLnStderr $ "Array index should be of type int but got " ++ actualType
  printCode code [pos]
printErrHolder code (InvalidArrayAccess pos notAnArrayType) = do
  putStrLnStderr $ "Error: Invalid array access at " ++ posToString pos
  putStrLnStderr $ "Expression is not an array but " ++ notAnArrayType
  printCode code [pos]
printErrHolder code (InvalidClassFieldAccess pos notAClassType) = do
  putStrLnStderr $ "Error: Invalid class field access at " ++ posToString pos
  putStrLnStderr $ "Expression is not a class but " ++ notAClassType
  printCode code [pos]
printErrHolder code (UnknownClassField fieldName className pos) = do
  putStrLnStderr $ "Error: Unknown class field " ++ fieldName ++ " in class " ++ className ++ " at " ++ posToString pos
  printCode code [pos]
printErrHolder code (InvalidClassMethodAccess pos notAClassType) = do
  putStrLnStderr $ "Error: Invalid class method access at " ++ posToString pos
  putStrLnStderr $ "Expression is not a class but " ++ notAClassType
  printCode code [pos]
printErrHolder code (UnknownClassMethod methodName className pos) = do
  putStrLnStderr $ "Error: Unknown class method " ++ methodName ++ " in class " ++ className ++ " at " ++ posToString pos
  printCode code [pos]
printErrHolder code (InvalidNegType pos actualType) = do
  putStrLnStderr $ "Error: Invalid negation type at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type int but got " ++ actualType
  printCode code [pos]
printErrHolder code (InvalidNotType pos actualType) = do
  putStrLnStderr $ "Error: Invalid negation type at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type bool but got " ++ actualType
  printCode code [pos]
printErrHolder code (InvalidMulType pos type1 type2) = do
  putStrLnStderr $ "Error: Invalid multiplication type at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type int but got " ++ type1 ++ " and " ++ type2
  printCode code [pos]
printErrHolder code (InvalidAddType pos type1 type2) = do
  putStrLnStderr $ "Error: Invalid addition type at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type int but got " ++ type1 ++ " and " ++ type2
  printCode code [pos]
printErrHolder code (InvalidRelType pos type1 type2) = do
  putStrLnStderr $ "Error: Invalid relation type at " ++ posToString pos
  putStrLnStderr $ "Expression should be of the same type but got " ++ type1 ++ " and " ++ type2
  printCode code [pos]
printErrHolder code (InvalidRelOp pos type1 type2) = do
  putStrLnStderr $ "Error: Invalid relation operator at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type int but got " ++ type1 ++ " and " ++ type2
  printCode code [pos]
printErrHolder code (InvalidRelTypeVoid pos type1 type2) = do
  putStrLnStderr $ "Error: Invalid relation type at " ++ posToString pos
  putStrLnStderr $ "Expression should not compare voids but got " ++ type1 ++ " and " ++ type2
  printCode code [pos]
printErrHolder code (InvalidAndType pos type1 type2) = do
  putStrLnStderr $ "Error: Invalid and type at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type bool but got " ++ type1 ++ " and " ++ type2
  printCode code [pos]
printErrHolder code (InvalidOrType pos type1 type2) = do
  putStrLnStderr $ "Error: Invalid or type at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type bool but got " ++ type1 ++ " and " ++ type2
  printCode code [pos]
printErrHolder code (CannotAssignType pos type1 type2) = do
  putStrLnStderr $ "Error: Cannot assign type at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type " ++ type1 ++ " but got " ++ type2
  printCode code [pos]
printErrHolder code (CannotIncrementNotInt pos type1) = do
  putStrLnStderr $ "Error: Cannot increment not int at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type int but got " ++ type1
  printCode code [pos]
printErrHolder code (CannotDecrementNotInt pos type1) = do
  putStrLnStderr $ "Error: Cannot decrement not int at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type int but got " ++ type1
  printCode code [pos]
printErrHolder code (IfCondShouldBeBoolean pos type1) = do
  putStrLnStderr $ "Error: If condition should be boolean at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type bool but got " ++ type1
  printCode code [pos]
printErrHolder code (WhileCondShouldBeBoolean pos type1) = do
  putStrLnStderr $ "Error: While condition should be boolean at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type bool but got " ++ type1
  printCode code [pos]
printErrHolder code (ForIterableNotArray pos type1) = do
  putStrLnStderr $ "Error: For iterable should be array at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type array but got " ++ type1
  printCode code [pos]
printErrHolder code (InvalidNullCast pos type1) = do
  putStrLnStderr $ "Error: Invalid null cast at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type class but got " ++ type1
  printCode code [pos]
printErrHolder code (CannotChangeNotLValue pos) = do
  putStrLnStderr $ "Error: Cannot change not lvalue at " ++ posToString pos
  putStrLnStderr "Expression should be lvalue"
  printCode code [pos]
printErrHolder code (IntOverflow pos value) = do
  putStrLnStderr $ "Error: Int overflow at " ++ posToString pos
  putStrLnStderr $ "Expression should be in range [-2^31, 2^31-1] but got " ++ show value
  printCode code [pos]
printErrHolder code (ForIterableTypeMismatch pos type1 type2) = do
  putStrLnStderr $ "Error: For iterable type mismatch at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type " ++ type1 ++ " but got " ++ type2
  printCode code [pos]
printErrHolder code (DeclTypeMismatch pos type1 type2) = do
  putStrLnStderr $ "Error: Declaration type mismatch at " ++ posToString pos
  putStrLnStderr $ "Expression should be of type " ++ type1 ++ " but got " ++ type2
  printCode code [pos]
printErrHolder code (InvalidFieldName pos fieldName) = do
  putStrLnStderr $ "Error: Invalid field name at " ++ posToString pos
  putStrLnStderr $ "Field name should be valid identifier but got " ++ fieldName
  printCode code [pos]
printErrHolder code (DuplicateFunctionName funName poses) = do
  putStrLnStderr $ "Error: Duplicate function name: " ++ funName
  mapM_ (putStrLnStderr . ("  " ++) . posToString) poses
  printCode code poses
printErrHolder code (DuplicateBuiltinFunctionName funName poses) = do
  putStrLnStderr $ "Error: Duplicate builtin function name: " ++ funName
  mapM_ (putStrLnStderr . ("  " ++) . posToString) poses
  printCode code poses
printErrHolder code (ReturnVoid pos) = do
  putStrLnStderr $ "Error: Return void at " ++ posToString pos
  putStrLnStderr "Function cannot return void (even if is declared as void)"
  printCode code [pos]
printErrHolder code (DivisionByZero pos) = do
  putStrLnStderr $ "Error: Division by zero at " ++ posToString pos
  putStrLnStderr "Expression should not be zero"
  printCode code [pos]
printErrHolder code (InvalidReturnType pos actualType expectedType) = do
  putStrLnStderr $ "Error: Invalid return type at " ++ posToString pos
  putStrLnStderr $ "Function should return a value of type " ++ expectedType ++ " but returns " ++ actualType
  printCode code [pos]
printErrHolder code (DuplicateClassName className poses) = do
  putStrLnStderr $ "Error: Duplicate class name: " ++ className
  mapM_ (putStrLnStderr . ("  " ++) . posToString) poses
  printCode code poses
printErrHolder code (InvalidSuperClassMethodRetType methodName currentPos currentRetType superPos superRetType) = do
  putStrLnStderr $ "Error: Invalid superclass method return type at " ++ posToString currentPos
  putStrLnStderr $ "Method " ++ methodName ++ " should return a value of type " ++ superRetType ++ " but returns " ++ currentRetType
  printCode code [currentPos, superPos]
printErrHolder code (InvalidSuperClassMethodArgs methodName currentPos currentArgs superPos superArgs) = do
  putStrLnStderr $ "Error: Invalid superclass method arguments at " ++ posToString currentPos
  putStrLnStderr $ "Method " ++ methodName ++ " should have arguments of types (" ++ intercalate ", " superArgs ++ ") but got (" ++ intercalate ", " currentArgs ++ ")"
  printCode code [currentPos, superPos]
printErrHolder code (InvalidClassName pos className) = do
  putStrLnStderr $ "Error: Invalid class name at " ++ posToString pos
  putStrLnStderr $ "Class name should be valid identifier but got " ++ className
  printCode code [pos]
printErrHolder code (InvalidArrayType pos arrayType) = do
  putStrLnStderr $ "Error: Invalid array type at " ++ posToString pos
  putStrLnStderr $ "Array type should be valid identifier but got " ++ arrayType
  printCode code [pos]

-- code printing ---------------------------------------------------------------

printCode :: String -> [Pos] -> IO ()
printCode = printCodeInternal 3

printCodeInternal :: Int -> String -> [Pos] -> IO ()
printCodeInternal offset sourceCode elementsToMark = do
  let codeLines = lines sourceCode
  let minCodeLine = 1
  let maxCodeLine = length codeLines
  let codeRanges = mergeCodeRanges $ makeCodeRanges minCodeLine maxCodeLine offset elementsToMark
  let codeWithLineNumbers = zip [minCodeLine ..] codeLines
  let codesInRange = takeOnlyInRange codeRanges codeWithLineNumbers
  let maxNumberLength = length $ show maxCodeLine
  let markers = makeMarkers $ normalizeMarkers elementsToMark
  let merged = sort $ codesInRange ++ markers
  putStrLnStderr ""
  mapM_ (printCodeLine maxNumberLength) merged
  where
    lineNumberToString :: Int -> Int -> String
    lineNumberToString lineNumber maxNumberLength =
      let lineNumberString = show lineNumber
       in let lineNumberStringLength = length lineNumberString
           in let numberOfSpaces = maxNumberLength - lineNumberStringLength
               in replicate numberOfSpaces ' ' ++ lineNumberString
    printCodeLine :: Int -> (Int, Int, String) -> IO ()
    printCodeLine maxNumberLength (lineNumber, 0, code) =
      putStrLnStderr $ lineNumberToString lineNumber maxNumberLength ++ ": " ++ code
    printCodeLine maxNumberLength (_, 1, code) =
      putStrLnStderr $ replicate maxNumberLength '-' ++ "- " ++ code
    printCodeLine maxNumberLength (_, _, code) =
      putStrLnStderr $ replicate maxNumberLength '=' ++ "= " ++ code

makeCodeRanges :: Int -> Int -> Int -> [Pos] -> [(Int, Int)]
makeCodeRanges minCodeLine maxCodeLine offset positions =
  [(max (line - offset) minCodeLine, min (line + offset) maxCodeLine) | Just (line, _) <- positions]

mergeCodeRanges :: [(Int, Int)] -> [(Int, Int)]
mergeCodeRanges x | length x <= 1 = x
mergeCodeRanges codeRanges =
  let x : xs = sort codeRanges
   in sort $ mergeCodeRangesInner [x] xs
  where
    mergeCodeRangesInner acc ranges =
      case ranges of
        [] -> acc
        (start, end) : xs ->
          let (lastStart, lastEnd) : acc' = acc
           in if start <= lastEnd + 1
                then mergeCodeRangesInner ((lastStart, max lastEnd end) : acc') xs
                else mergeCodeRangesInner ((start, end) : acc) xs

takeOnlyInRange :: [(Int, Int)] -> [(Int, String)] -> [(Int, Int, String)]
takeOnlyInRange codeRanges codeWithLineNumbers =
  let codeRanges' = sort codeRanges
   in let codeWithLineNumbers' = sort codeWithLineNumbers
       in takeOnlyInRangeInner codeRanges' codeWithLineNumbers'
  where
    prepend :: Int -> [(Int, String)] -> [(Int, Int, String)]
    prepend w = map (\(lineNumber, code) -> (lineNumber, w, code))
    takeOnlyInRangeInner codeRangesSorted codeWithLineNumbersSorted =
      case codeRangesSorted of
        [] -> []
        (start, end) : xs ->
          let codeWithLineNumbers' = dropWhile (\(lineNumber, _) -> lineNumber < start) codeWithLineNumbersSorted
           in let codeWithLineNumbers'' = takeWhile (\(lineNumber, _) -> lineNumber <= end) codeWithLineNumbers'
               in prepend 0 codeWithLineNumbers'' ++ [(end, 2, "")] ++ takeOnlyInRangeInner xs codeWithLineNumbers'

makeMarkers :: [(Int, [Int])] -> [(Int, Int, String)]
makeMarkers = map (\(l, c) -> (l, 1, makeMarker c))
  where
    makeMarker :: [Int] -> String
    makeMarker cols =
      let cols' = sort $ nub cols
       in fst $ foldl' go ("", 0) cols'
      where
        go :: (String, Int) -> Int -> (String, Int)
        go (acc, lastCol) col =
          let spaces = replicate (col - lastCol - 1) ' '
           in (acc ++ spaces ++ "^", col)

groupByKey :: [(Int, Int)] -> [(Int, [Int])]
groupByKey = map (\x -> (fst $ head x, map snd x)) . groupBy (\x y -> fst x == fst y)

normalizeMarkers :: [Pos] -> [(Int, [Int])]
normalizeMarkers positions = groupByKey [(line, col) | Just (line, col) <- positions]
