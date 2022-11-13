#include <cassert>
#include <iostream>
#include "fuzzy.h"

using std::cout;
using std::endl;

int main() {
  TriFuzzyNum num1(1, 2, 3);
  TriFuzzyNum num2(0.5, 0.25, 0.75);
  constexpr TriFuzzyNum num3(1, 1, 1);

  cout << num1 << endl;  // (1, 2, 3)
  cout << num2 << endl;  // (0.25, 0.5, 0.75)
  cout << num3 << endl;  // (1, 1, 1)

  assert(num1 + num2 == TriFuzzyNum(1.25, 2.5, 3.75));
  assert(num1 - num2 == TriFuzzyNum(0.25, 1.5, 2.75));
  assert(num1 * num2 == TriFuzzyNum(0.25, 1, 2.25));

  assert(num1 > num2);
  assert(num1 >= num2);
  assert(num1 != num2);
  assert(num2 < num1);
  assert(num2 <= num1);
  assert(num3 == crisp_number(1));
  static_assert(num3 == crisp_number(1));

  assert(num1 <=> num2 > 0);
  assert(num2 <=> num1 < 0);

  assert((num1 += num3) == TriFuzzyNum(2, 3, 4));
  assert((num1 -= num3) == TriFuzzyNum(1, 2, 3));
  assert((num1 *= num3) == num1);
  assert((num1 += crisp_zero) == num1);

  static_assert(num3.lower_value() == 1);
  static_assert(num3.modal_value() == 1);
  static_assert(num3.upper_value() == 1);

  TriFuzzyNumSet fn_set1({TriFuzzyNum(1, 2, 4), TriFuzzyNum(1, 3, 6)});
  assert(fn_set1.arithmetic_mean() == TriFuzzyNum(1, 2.5, 5));

  fn_set1.insert(TriFuzzyNum(1, 4, 5));
  assert(fn_set1.arithmetic_mean() == TriFuzzyNum(1, 3, 5));

  fn_set1.remove(TriFuzzyNum(1, 3, 6));
  assert(fn_set1.arithmetic_mean() == TriFuzzyNum(1, 3, 4.5));

  TriFuzzyNumSet fn_set2({TriFuzzyNum(1, 2, 3),
                          TriFuzzyNum(1, 2, 3),
                          TriFuzzyNum(1, 5, 6)});
  assert(fn_set2.arithmetic_mean() == TriFuzzyNum(1, 3, 4));

  constinit static TriFuzzyNum num4(1.25, 2.25, 3.25);
  num4 += crisp_number(0.25);
  assert(num4 == TriFuzzyNum(1.5, 2.5, 3.5));
}
