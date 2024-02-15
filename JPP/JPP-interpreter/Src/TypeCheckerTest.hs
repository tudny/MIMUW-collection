import Src.Errors
import Src.Utils
import Src.Types
import Src.Jabba.Par
import Src.Jabba.Abs
import System.Directory.Internal.Prelude (exitFailure)
import Src.TypeChecker
import qualified Data.Map as Map

type TestCase = (String, String, Env, ErrType)

type FileTestCase = (String, ErrType)
type InlineTestCase = (String, Env, ErrType)

testFiles :: [FileTestCase]
testFiles = [
    ("12-000-zero", ZeroLiternalDiv),
    ("12-001-const", ImmutVar (Ident "x")),
    ("12-002-multiname", VarAlreadyDecl (Ident "x")),
    ("12-003-wrong-type", WrongType (Ident "s") VTString [VTInt]),
    ("12-004-wrong-type-var", WrongType (Ident "s") VTString [VTInt]),
    ("12-005-multi-wrong-type", WrongType (Ident "x") VTInt [VTString]),
    ("12-006-multi-wrong-type-var", WrongType (Ident "x") VTInt [VTString]),
    ("12-007-undef", NotDeclVar (Ident "x")),
    ("12-008-undef-op", NotDeclVar (Ident "x")),
    ("12-009-boolen", WrongType (Ident "x") VTInt [VTBool]),
    ("12-010-boolen2", WrongType (Ident "x") VTString [VTBool]),
    ("12-011-top-level-return", TopLevelProgramReturn VTInt),
    ("12-012-top-level-maybe-return", TopLevelProgramMaybeReturn VTInt),
    ("12-013-top-level-return2", TopLevelProgramReturn VTInt),
    ("12-014-top-level-maybe-return2", TopLevelProgramMaybeReturn VTInt),
    ("12-015-top-level-return3", TopLevelProgramReturn VTInt),
    ("12-016-top-level-return4", TopLevelProgramReturn VTInt),
    ("12-017-top-level-loop-flow", TopLevelProgramLoopFlow),
    ("12-018-top-level-loop-flow2", TopLevelProgramLoopFlow),
    ("12-019-top-level-loop-flow3", NotDeclVar (Ident "x")),
    ("12-020-top-level-loop-flow4", NotDeclVar (Ident "x")),
    ("12-021-top-level-loop-flow5", NotDeclVar (Ident "x")),
    ("12-022-top-level-return-flow", TopLevelProgramReturn VTInt),
    ("12-023-fun-loop-flow", FunctionLoopFlow (Ident "foo")),
    ("12-024-fun-loop-flow2", FunctionLoopFlow (Ident "foo")),
    ("12-025-fun-var-vis", NotDeclVar (Ident "x")),
    ("12-026-fun-var-vis2", NotDeclVar (Ident "x")),
    ("12-027-fun-var-vis3", NotDeclVar (Ident "x")),
    ("12-028-fun-return-match", MismatchedReturnTypes VTInt VTString),
    ("12-029-fun-return-match2", FunctionReturnMismatch (Ident "foo") VTInt VTString),
    ("12-030-fun-call-mut", ImmutMutPass 0),
    ("12-031-lambda-mismatch", WrongType (Ident "f") (Fn [(VTInt, VMConst, VRRef)] VTInt) [Fn [(VTInt, VMMut, VRRef)] VTInt]),
    ("12-032-function-hide", ConstantAssign (Ident "foo")),
    ("12-033-function-recursion", NotDeclVar (Ident "not_decl")),
    ("12-034-lambda-passing", NotDeclVar (Ident "not_decl")),
    ("12-035-const-in-fun", ConstantAssign (Ident "x")),
    ("12-036-function-wo-body", FunctionWithoutInitializer (Ident "f"))
  ]

inlineTests :: [InlineTestCase]
inlineTests = [
    (
      "x++;", 
      Env (Map.fromList [(Ident "x", (VTInt, VMConst))]), 
      ImmutVar (Ident "x")
    ),
    (
      "f(1);",
      Env (Map.fromList [(Ident "f", (Fn [] VTInt, fMut))]),
      WrongNumberOfArgs 0 1
    ),
    (
      "f();",
      Env (Map.fromList [(Ident "f", (Fn [(VTInt, VMConst, VRRef)] VTInt, fMut))]),
      WrongNumberOfArgs 1 0
    ),
    (
      "f(1, 2, 3, 4, 5, 6);",
      Env (Map.fromList [(Ident "f", (Fn [(VTInt, VMConst, VRRef)] VTInt, fMut))]),
      WrongNumberOfArgs 1 6
    ),
    (
      "f(\"abc\");",
      Env (Map.fromList [(Ident "f", (Fn [(VTInt, VMConst, VRRef)] VTInt, fMut))]),
      WrongTypeArg 0 VTString VTInt
    ),
    (
      "f(1, \"abc\");",
      Env (Map.fromList [(Ident "f", (Fn [(VTInt, VMConst, VRCopy), (VTInt, VMConst, VRCopy)] VTInt, fMut))]),
      WrongTypeArg 1 VTString VTInt
    ),
    (
      "f(x);",
      Env (Map.fromList [(Ident "f", (Fn [(VTInt, VMConst, VRRef)] VTInt, fMut))]),
      NotDeclVar (Ident "x")
    ),
    (
      "f(10);",
      Env (Map.fromList [(Ident "f", (Fn [(VTInt, VMMut, VRRef)] VTInt, fMut))]),
      ExprRefPass 0
    ),
    (
      "f(x);",
      Env (
        Map.fromList [(Ident "x", (VTInt, VMConst)), (Ident "f", (Fn [(VTInt, VMMut, VRRef)] VTInt, fMut))]
      ),
      ImmutMutPass 0
    ),
    (
      "f(x + 1);",
      Env (
        Map.fromList [(Ident "x", (VTInt, VMConst)), (Ident "f", (Fn [(VTInt, VMMut, VRRef)] VTInt, fMut))]
      ),
      ExprRefPass 0
    ),
    (
      "f(x + 0);",
      Env (
        Map.fromList [(Ident "x", (VTInt, VMConst)), (Ident "f", (Fn [(VTInt, VMMut, VRRef)] VTInt, fMut))]
      ),
      ExprRefPass 0
    ),
    (
      "-\"abc\";",
      Env Map.empty,
      WrongTypeOp "negation" VTString
    ),
    (
      "!\"abc\";",
      Env Map.empty,
      WrongTypeOp "negation" VTString
    ),
    (
      "-true;",
      Env Map.empty,
      WrongTypeOp "negation" VTBool
    ),
    (
      "!1;",
      Env Map.empty,
      WrongTypeOp "negation" VTInt
    ),
    (
      "1+\"abc\";",
      Env Map.empty,
      WrongTypeBiOp "addition" VTInt VTString
    ),
    (
      "1-\"abc\";",
      Env Map.empty,
      WrongTypeBiOp "subtraction" VTInt VTString
    ),
    (
      "\"abc\"+1;",
      Env Map.empty,
      WrongTypeBiOp "addition" VTString VTInt
    ),
    (
      "\"abc\"-1;",
      Env Map.empty,
      WrongTypeBiOp "subtraction" VTString VTInt
    ),
    (
      "\"abc\"*1;",
      Env Map.empty,
      WrongTypeBiOp "multiplication" VTString VTInt
    ),
    (
      "\"abc\"/1;",
      Env Map.empty,
      WrongTypeBiOp "division" VTString VTInt
    ),
    (
      "\"abc\"%1;",
      Env Map.empty,
      WrongTypeBiOp "modulo" VTString VTInt
    ),
    (
      "\"abc\"<1;",
      Env Map.empty,
      WrongTypeBiOp "less than" VTString VTInt
    ),
    (
      "\"abc\">1;",
      Env Map.empty,
      WrongTypeBiOp "greater than" VTString VTInt
    ),
    (
      "\"abc\"<=1;",
      Env Map.empty,
      WrongTypeBiOp "less or equal" VTString VTInt
    ),
    (
      "\"abc\">=1;",
      Env Map.empty,
      WrongTypeBiOp "greater or equal" VTString VTInt
    ),
    (
      "\"abc\"==1;",
      Env Map.empty,
      WrongTypeBiOp "equality" VTString VTInt
    ),
    (
      "\"abc\"!=1;",
      Env Map.empty,
      WrongTypeBiOp "inequality" VTString VTInt
    ),
    (
      "\"abc\"&&1;",
      Env Map.empty,
      WrongTypeBiOp "and" VTString VTInt
    ),
    (
      "\"abc\"||1;",
      Env Map.empty,
      WrongTypeBiOp "or" VTString VTInt
    ),
    (
      "1+true;",
      Env Map.empty,
      WrongTypeBiOp "addition" VTInt VTBool
    ),
    (
      "1-true;",
      Env Map.empty,
      WrongTypeBiOp "subtraction" VTInt VTBool
    ),
    (
      "true+1;",
      Env Map.empty,
      WrongTypeBiOp "addition" VTBool VTInt
    ),
    (
      "true-1;",
      Env Map.empty,
      WrongTypeBiOp "subtraction" VTBool VTInt
    ),
    (
      "true*1;",
      Env Map.empty,
      WrongTypeBiOp "multiplication" VTBool VTInt
    ),
    (
      "true/1;",
      Env Map.empty,
      WrongTypeBiOp "division" VTBool VTInt
    ),
    (
      "true%1;",
      Env Map.empty,
      WrongTypeBiOp "modulo" VTBool VTInt
    ),
    (
      "true<1;",
      Env Map.empty,
      WrongTypeBiOp "less than" VTBool VTInt
    ),
    (
      "true>1;",
      Env Map.empty,
      WrongTypeBiOp "greater than" VTBool VTInt
    ),
    (
      "true<=1;",
      Env Map.empty,
      WrongTypeBiOp "less or equal" VTBool VTInt
    ),
    (
      "true>=1;",
      Env Map.empty,
      WrongTypeBiOp "greater or equal" VTBool VTInt
    ),
    (
      "true==1;",
      Env Map.empty,
      WrongTypeBiOp "equality" VTBool VTInt
    ),
    (
      "true!=1;",
      Env Map.empty,
      WrongTypeBiOp "inequality" VTBool VTInt
    ),
    (
      "true?1:\"abc\";",
      Env Map.empty,
      TernaryMismatch VTInt VTString
    ),
    (
      "true?1:\"abc\";",
      Env Map.empty,
      TernaryMismatch VTInt VTString
    ),
    (
      "var x: Boolean = (a < b); z;",
      Env (Map.fromList [
        (Ident "a", (VTInt, VMConst)),
        (Ident "b", (VTInt, VMConst))
        ]),
      NotDeclVar (Ident "z")
    ),
    (
      "var x: Boolean = (!b); z;",
      Env (Map.fromList [
        (Ident "b", (VTBool, VMConst))
        ]),
      NotDeclVar (Ident "z")
    ),
    (
      "(x && y || z) ? 1 : 2;",
      Env (Map.fromList [
        (Ident "x", (VTBool, VMConst)),
        (Ident "y", (VTBool, VMConst)),
        (Ident "z", (VTInt, VMConst))
        ]),
      WrongTypeBiOp "or" VTBool VTInt
    ),
    (
      "(x && y) ? \"ERROR\" : 1;",
      Env (Map.fromList [
        (Ident "x", (VTBool, VMConst)),
        (Ident "y", (VTBool, VMConst))
        ]),
      TernaryMismatch VTString VTInt
    ),
    (
      "1 ? 2 : 3;",
      Env Map.empty,
      WrongTypeOp "ternary operator condition" VTInt
    ),
    (
      "x = 1;",
      Env (Map.fromList [
        (Ident "x", (VTInt, VMConst))
        ]),
      ConstantAssign (Ident "x")
    ),
    (
      "x = 1;",
      Env (Map.fromList [
        (Ident "x", (VTString, VMMut))
        ]),
      WrongType (Ident "x") VTString [VTInt]
    ),
    (
      "val x: Integer = 1; x = 2;",
      Env Map.empty,
      ConstantAssign (Ident "x")
    ),
    (
      "var x: String = \"abc\"; x = 2;",
      Env Map.empty,
      WrongType (Ident "x") VTString [VTInt]
    ),
    (
      "var x: Boolean = true; x = 2;",
      Env Map.empty,
      WrongType (Ident "x") VTBool [VTInt]
    ),
    (
      "var x: Integer = 1; x = \"abc\";",
      Env Map.empty,
      WrongType (Ident "x") VTInt [VTString]
    ),
    (
      "var x: Integer = 1; x = true;",
      Env Map.empty,
      WrongType (Ident "x") VTInt [VTBool]
    ),
    (
      "var b = || { return true; }; b = || { return 1; };",
      Env Map.empty,
      WrongType (Ident "b") (Fn [] VTBool) [Fn [] VTInt]
    )
  ]


checkErrorMatch :: ErrType -> ErrHolder -> Bool
checkErrorMatch t h = 
    case h of
        TypeChecker _ t' -> t == t'
        _ -> False


testSingleInlineCase :: TestCase -> IO ()
testSingleInlineCase (name, content, env, expected) = do
  let ts = myLexer content
  case go ts of
    Left err -> do
      if checkErrorMatch expected err
        then putStrLn $ "Success " ++ "\"" ++ trim name ++ "\""
        else do
          putStrLn $ "Failure " ++ show name
          putStrLn $ "Expected: " ++ show expected
          putStrLn $ "Got: " ++ show err
          exitFailure
    Right _ -> do
      putStrLn $ "Failure " ++ show name
      putStrLn $ "Expected: " ++ show expected
      putStrLn $ "Got: " ++ "No error"
      exitFailure
  where
    trim x | length x > 40 = take 40 x ++ "..."
           | otherwise = x
    go ts = do
      tree <- left ParserErr $ pProgram ts
      typeCheckWithEnv env tree
      pure ()


fileCaseResolve :: FileTestCase -> IO TestCase
fileCaseResolve (name, expected) = do
  fileContent <- readFile ("bad/" ++ name ++ ".jbb")
  pure (name, fileContent, emptyEnv, expected)

inlineCaseResolve :: InlineTestCase -> IO TestCase
inlineCaseResolve (content, env, expected) = pure (content, content, env, expected)


main :: IO ()
main = do
  putStrLn "Running tests..."
  putStrLn $ "Detected files number: " ++ show (length testFiles)
  putStrLn $ "Detected inline tests number: " ++ show (length inlineTests) 
  fileCases <- mapM fileCaseResolve testFiles
  inlineCases <- mapM inlineCaseResolve inlineTests
  mapM_ testSingleInlineCase $ fileCases ++ inlineCases
