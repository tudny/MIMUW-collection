-- Zadanie 2 Poniżej zaproponowana jest składnia pewnego języka programowania

type Var = String
type PName = String
type Number = Int
data Expr = Number Int | Variable Var | Plus Expr Expr | Minus Expr Expr | Mul Expr Expr
data Instr = Skip | Assign Var Expr | Seq Instr Instr | If Expr Instr Instr | ProcDef PName Var Instr | ProcCall PName Expr

type Loc = Int
type Store = Loc -> Number
type EnvV = Var -> Loc
type EnvP = PName -> Proc
type Proc = Expr -> EnvV -> (EnvP, Store) -> (EnvP, Store)

-- Wyrażenia i instrukcje oprócz wymienionych poniżej mają standardowe znaczenie. Można
-- założyć, że wszystkie zmienne indywiduowe x, y, . . . mają początkowo określone wartości (są
-- zainicjowane). Wiązanie zmiennych indywiduowych jest statyczne, dopuszczalne wartości tych
-- zmiennych są typu Num. Wyrażenia arytmetyczne wyliczają się do wartości liczbowych Num
-- bez efektów ubocznych. Można założyć, że dana jest semantyka dla wyrażeń w stylu denotacyjnym, należy jedynie jawnie zadeklarować (sensowny) typ funkcji semantycznej. Podobnie można przyjąć, że dostępny jest nieskończony zbiór lokacji Loc oraz funkcja matematyczna newloc: (Loc *fin Num) → Loc spełniająca newloc(s) ∈/ dom(s) dla każdego
-- s: Loc *fin Num (gdzie Loc *fin Num to zbiór funkcji częściowych z Loc do Num, których dziedzina jest skończona).
-- Język jest wzbogacony o jednoparametrowe procedury. Definicja procedury następuje z użyciem instrukcji P := proc(x){I}, która przypisuje nazwie procedury P procedurę o ciele I z
-- parametrem formalnym x. Zmienne indywiduowe występujące w ciele I tak zdefiniowanej procedury powinny zostać w tym momencie związane statycznie, niezależnie od kontekstu, w jakim
-- procedura będzie później wywoływana. Zmienna x będąca parametrem formalnym procedury
-- staje się zmienną lokalną widoczną w ciele I (i w procedurach w nim wprowadzonych). Zmienna
-- ta przesłania wcześniej zdefiniowaną zmienną o tej samej nazwie.
-- Wywołanie tak zdefiniowanej procedury P instrukcją call P(e) powoduje obliczenie aktualnej
-- wartości wyrażenia e, a następnie przypisanie tej wartości zmiennej x wprowadzonej dla tego
-- wywołania procedury i wykonanie ciała I. Oznacza to, że przekazywanie parametrów odbywa
-- się przez wartość.
-- Instrukcja przypisania procedur P := P1 then P2 zmienia globalną wartość procedury P
-- w ten sposób, że staje się ona złożeniem procedur P1 i P2 (identyfikatory procedur P, P1 i
-- P2 nie muszą być wzajemnie różne). Wywołanie tak powstałej procedury instrukcją call P(e)
-- powoduje wywołanie po kolei procedur P1(e) i P2(e), zgodnie z ich wartościami sprzed przypisania. W szczególności wyrażenie e będzie obliczane (przynajmniej) dwukrotnie, raz na potrzeby
-- wykonania procedury P1, a następnie na potrzeby wykonania P2.
-- Wiązanie identyfikatorów procedur jest dynamiczne; w szczególności przy wykonaniu ciała
-- procedury może dojść do rekurencyjnego wywołania tejże procedury lub innej procedury, która
-- ją wywoła. Początkowo wszystkie nazwy procedur są zainicjowane na procedury, które nic nie
-- robią, czyli P := proc(x){skip}.
-- Zadanie
-- Zadanie polega na napisaniu semantyki denotacyjnej dla kategorii syntaktycznej Inst instrukcji
-- powyższego języka. W tym celu należy jednoznacznie zdefiniować typy funkcji semantycznych
-- dla Inst oraz Expr wraz ze wszystkimi używanymi typami pomocniczymi. Należy też podać
-- równania semantyczne dla wszystkich postaci instrukcji Inst w tym języku.
-- verte!

-- Uwaga
-- Przypisywanie identyfikatorów procedur w tym języku odbywa się w sposób dynamiczny. Powoduje to, że faktyczna wartość procedury zmienia się w dynamicznie, w trakcie wykonywania
-- programu. W szczególności, każde wywołanie procedury o nazwie P powoduje jej wywołanie
-- zgodnie z aktualną wartością tego identyfikatora. W odróżnieniu od zmiennych indywiduowych,
-- które mają zakresy widoczności i globalny, i lokalne, zmienne procedurowe mają zawsze zakres
-- widoczności globalny, co oznacza, że zmiana wartości identyfikatora procedury w ciele jakiejś
-- procedury ma efekt również po wyjściu z tej ostatniej. Jest to niestandardowa cecha tego języka
-- — w przypadku większości języków rozważanych na tym przedmiocie wiązanie procedur odbywa
-- się w sposób statyczny, a ich wartości są niezmienne w czasie.


newloc :: EnvV -> Loc
newloc s = 1 + maximum (map fst (filter (\(x,y) -> y /= 0) (assocs s)))

fE :: Expr -> EnvV -> Store -> Number
fI :: Instr -> EnvV -> (EnvP, Store) -> (EnvP, Store)


fE (Number n) _ _ = n
fI Skip _ r = r


