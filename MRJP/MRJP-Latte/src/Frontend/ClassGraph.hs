module Frontend.ClassGraph where

import Control.Monad.State
import Data.Graph
import qualified Data.Map as Map
import Data.Maybe
import qualified Data.Set as Set
import Data.Tuple.Extra
import Data.Tuple.Lazy
import Debug.Trace (trace)
import Frontend.AbstractSyntaxTree

sortClasses :: [(String, Maybe String)] -> [String]
sortClasses allClasses =
  let edges = map (\(x, y) -> (x, x, maybeToList y)) allClasses
   in let (graph, vertexToNode, _) = graphFromEdges edges
       in let sorted = topSort graph
           in reverse $ map (fst3 . vertexToNode) sorted

enrichClassesWithSuperClasses :: Map.Map String ASTClass -> Map.Map String ASTClass
enrichClassesWithSuperClasses classes =
  let sortedClassNames = sortClasses $ toNodes classes
   in foldl enrich Map.empty sortedClassNames
  where
    toNodes :: Map.Map String ASTClass -> [(String, Maybe String)]
    toNodes classes =
      map
        ( \(name, ASTClassObj {className = className, classSuper = classSuper}) -> (className, classSuper)
        )
        $ Map.toList classes
    enrich :: Map.Map String ASTClass -> String -> Map.Map String ASTClass
    enrich enriched currentClass =
      let ASTClassObj {classSuper = classSuper} = classes Map.! currentClass
       in case classSuper of
            Nothing -> Map.insert currentClass (classes Map.! currentClass) enriched
            Just superClass ->
              let ASTClassObj {classFields = classFields, classMethods = classMethods} = classes Map.! currentClass
                  ASTClassObj {classFields = superFields, classMethods = superMethods} = enriched Map.! superClass
                  mergedFields = superFields ++ classFields
                  mergedMethods = mergeMethods currentClass superMethods classMethods
               in Map.insert
                    currentClass
                    ( ASTClassObj
                        { className = currentClass,
                          classFields = mergedFields,
                          classMethods = mergedMethods,
                          classSuper = classSuper
                        }
                    )
                    enriched

mergeMethods :: String -> [ASTMethod] -> [ASTMethod] -> [ASTMethod]
mergeMethods className parentMethods myMethods = do
  let myFunsAsSet = Set.fromList $ map methodFun myMethods
      (parentMethodsEnriched, myOwnMethods) = runState (mapM (mergeMethod className) parentMethods) myFunsAsSet
      myOwnMethodsEnriched = map (\fun -> ASTMethodObj {methodFun = fun, methodIImplement = True, methodImplementor = className}) (Set.elems myOwnMethods)
  parentMethodsEnriched ++ myOwnMethodsEnriched
  where
    mergeMethod :: String -> ASTMethod -> State (Set.Set ASTFun) ASTMethod
    mergeMethod className me@ASTMethodObj {methodFun = funPart, methodImplementor = implementor} = do
      doIImplementIt <- gets $ Set.member funPart
      let actualImplementor = if doIImplementIt then className else implementor
      when doIImplementIt $ modify $ Set.delete funPart
      pure $ ASTMethodObj {methodFun = funPart, methodIImplement = doIImplementIt, methodImplementor = actualImplementor}

-- >>> let classA = ASTClassObj {className = "A", classFields = Map.fromList [("a", ASTFieldObj {classFieldName = "a", classFieldType = ASTInt})], classMethods = Map.fromList [("a", ASTFunObj {funName = "a", funArgs = [], funRet = ASTInt})], classSuper = Nothing}
-- >>> let classB = ASTClassObj {className = "B", classFields = Map.fromList [("b", ASTFieldObj {classFieldName = "b", classFieldType = ASTInt})], classMethods = Map.fromList [("b", ASTFunObj {funName = "b", funArgs = [], funRet = ASTInt})], classSuper = Just "A"}
-- >>> let classes = Map.fromList [("A", classA), ("B", classB)]
-- >>> enrichClassesWithSuperClasses classes
-- fromList [("A",ASTClassObj {className = "A", classFields = fromList [("a",ASTFieldObj {classFieldName = "a", classFieldType = ASTInt})], classMethods = fromList [("a",ASTFunObj {funName = "a", funArgs = [], funRet = ASTInt})], classSuper = Nothing}),("B",ASTClassObj {className = "B", classFields = fromList [("a",ASTFieldObj {classFieldName = "a", classFieldType = ASTInt}),("b",ASTFieldObj {classFieldName = "b", classFieldType = ASTInt})], classMethods = fromList [("a",ASTFunObj {funName = "a", funArgs = [], funRet = ASTInt}),("b",ASTFunObj {funName = "b", funArgs = [], funRet = ASTInt})], classSuper = Just "A"})]
