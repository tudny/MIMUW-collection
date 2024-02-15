
import System.Environment
import System.IO (withFile, IOMode (ReadMode), hGetContents)

-- a. Napisz program który wypisze swoje argumenty, każdy w osobnej linii (wskazówka: System.Environment.getArgs)

main1 = do
    progName <- getProgName
    args <- getArgs
    putStrLn $ "Running: " ++ progName
    putStrLn (unlines args)

-- b. Napisz program, który będzie pytał użytkownika o ulubiony język programowania, tak długo aż odpowiedzią będzie 'Haskell' ;)

main2 = do
    putStrLn "What is your favorite programming language?"
    lang <- getLine
    if lang == "Haskell"
        then putStrLn "Great!"
        else main2

-- c. Napisz uproszczoną wersję programu wc (wypisującą liczbę linii, słów i znaków w pliku o nazwie podanej jako argument, bądź stdin jeśli bez argumentu).

main = do
    args <- getArgs
    let fileName = if null args then "" else head args
    withFile fileName ReadMode $ \handle -> do
        contents <- hGetContents handle
        let linesCount = length $ lines contents
        let wordsCount = length $ words contents
        let charsCount = length contents
        putStrLn $ "Lines: " ++ show linesCount
        putStrLn $ "Words: " ++ show wordsCount
        putStrLn $ "Chars: " ++ show charsCount

-- Wskazówki:

-- getLine :: IO String
-- putStrLn :: String -> IO ()

-- import System.Environment
-- getArgs :: IO [String]
-- getProgName :: IO String

-- import System.IO
-- withFile :: FilePath -> IOMode -> (Handle -> IO r) -> IO r
-- ReadMode :: IOMode
-- hGetContents :: Handle -> IO String



