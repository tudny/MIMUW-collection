#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>
#include <complex>
#include <algorithm>
#include <omp.h>

#include "utils/bmp.cpp"


void compress(const uint32_t valuesCount, const int accuracy,
    const uint8_t *values, float *Xreal, float *Ximag) {
  // values, Xreal and Ximag are values describing single color of single row of bitmap.
  // This function will be called once per each (color, row) combination.
  for (int k = 0; k < accuracy; k++) {
      for (int i = 0; i < valuesCount; i++) {
          float theta = (2 * M_PI * k * i) / valuesCount;
          Xreal[k] += values[i] * cos(theta);
          Ximag[k] -= values[i] * sin(theta);
      }
  }
}

void decompress(const uint32_t valuesCount, const int accuracy,
    uint8_t *values, const float *Xreal, const float *Ximag) {
  // values, Xreal and Ximag are values describing single color of single row of bitmap.
  // This function will be called once per each (color, row) combination.
  std::vector<float> rawValues(valuesCount, 0);

  for (int i = 0; i < valuesCount; i++) {
      for (int k = 0; k < accuracy; k++) {
          float theta = (2 * M_PI * k * i) / valuesCount;
          rawValues[i] += Xreal[k] * cos(theta) + Ximag[k] * sin(theta);
      }
      values[i] = rawValues[i] / valuesCount;
  }
}

void compressPar(const uint32_t valuesCount, const int accuracy,
    const uint8_t *values, float *Xreal, float *Ximag) {
  #pragma omp parallel for reduction(+:Xreal[:accuracy], Ximag[:accuracy]) collapse(2)
  for (int k = 0; k < accuracy; k++) {
      for (int i = 0; i < valuesCount; i++) {
          float theta = (2 * M_PI * k * i) / valuesCount;
          Xreal[k] += values[i] * cos(theta);
          Ximag[k] -= values[i] * sin(theta);
      }
  }
}

void decompressPar(const uint32_t valuesCount, const int accuracy,
    uint8_t *values, const float *Xreal, const float *Ximag) {

  std::vector<float> rawValues(valuesCount, 0);
  std::vector<std::vector<float>> cos_values(2, std::vector<float>(valuesCount, 1));
  std::vector<std::vector<float>> sin_values(2, std::vector<float>(valuesCount, 0));

  #pragma omp parallel
  {
    for (int k = 0; k < accuracy; k++) {
      int k_mod_2 = k % 2;
      #pragma omp master
      for (int i = 0; i < valuesCount; i++) {
        rawValues[i] += Xreal[k] * cos_values[k_mod_2][i] + Ximag[k] * sin_values[k_mod_2][i];
      }

      if (k + 1 < accuracy) {
        int k_plus_1_mod_2 = (k + 1) % 2;
        float m_pi_2_k_plus_1 = (2 * M_PI * (k + 1)) / valuesCount;
        #pragma omp for
        for (int i = 0; i < valuesCount; i++) {
          float theta = m_pi_2_k_plus_1 * i;
          cos_values[k_plus_1_mod_2][i] = cos(theta);
          sin_values[k_plus_1_mod_2][i] = sin(theta);
        }
      }

      #pragma omp barrier
    }

    #pragma omp master
    for (int i = 0; i < valuesCount; i++) {
      values[i] = rawValues[i] / valuesCount;
    }
  }
}

void compressNaivePar(const uint32_t valuesCount, const int accuracy,
    const uint8_t *values, float *Xreal, float *Ximag) {
  #pragma omp parallel for collapse(2)
  for (int k = 0; k < accuracy; k++) {
      for (int i = 0; i < valuesCount; i++) {
          float theta = (2 * M_PI * k * i) / valuesCount;
          #pragma omp atomic update
          Xreal[k] += values[i] * cos(theta);
          #pragma omp atomic update
          Ximag[k] -= values[i] * sin(theta);
      }
  }
}

int read_env_var(const char *name, int default_value) {
  char *env = getenv(name);
  if (env == NULL) {
    return default_value;
  }
  return atoi(env);
}

std::string make_filename(
  const std::string &prefix,
  const std::string &suffix,
  const size_t accuracy,
  const size_t threads
) {
  return "imgs/" + prefix + "_" + std::to_string(accuracy) + "_" + std::to_string(threads) + suffix;
}

int main() {
  size_t accuracy = read_env_var("ACCURACY", 16);
  size_t threads = read_env_var("OMP_NUM_THREADS", 1);
  BMP bmp;

  bmp.read("example.bmp"); // We are interested in values from range [8; 32]
  float seq_compress_time = bmp.compress(compress, accuracy);
  float seq_decompress_time = bmp.decompress(decompress);
  bmp.write(make_filename("example_result_seq", ".bmp", accuracy, threads));

  bmp.read("example.bmp");
  float par_fast_compress_time = bmp.compress(compressPar, accuracy);
  float __not__used__ = bmp.decompress(decompressPar);
  bmp.write(make_filename("example_result_par_fast", ".bmp", accuracy, threads));

  bmp.read("example.bmp");
  float par_naive_compress_time = bmp.compress(compressNaivePar, accuracy);
  float par_decompress_time = bmp.decompress(decompressPar);
  bmp.write(make_filename("example_result_par_naive", ".bmp", accuracy, threads));

  printf("SEQ_COMPRESS_TIME:%.4f\n", seq_compress_time);
  printf("SEQ_DECOMPRESS_TIME:%.4f\n", seq_decompress_time);
  printf("PAR_FAST_COMPRESS_TIME:%.4f\n", par_fast_compress_time);
  printf("PAR_NAIVE_COMPRESS_TIME:%.4f\n", par_naive_compress_time);
  printf("PAR_DECOMPRESS_TIME:%.4f\n", par_decompress_time);

  return 0;
}

