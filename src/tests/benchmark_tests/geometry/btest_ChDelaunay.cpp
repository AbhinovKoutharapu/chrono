// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================

#include <benchmark/benchmark.h>
#include <vector>
#include <random>

#include "chrono/geometry/ChDelaunay.h"
#include "chrono/geometry/ChDelaunay.cpp"

using namespace chrono;

// Helper to generate a cloud of random points
static std::vector<ChVector3d> GenerateRandomPoints(int n) {
    std::vector<ChVector3d> points;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(-10.0, 10.0);
    for (int i = 0; i < n; ++i) {
        points.emplace_back(dis(gen), dis(gen), dis(gen));
    }
    return points;
}

// Benchmark function
static void BM_Delaunay3D(benchmark::State& state) {
    auto points = GenerateRandomPoints((int)state.range(0));

    for (auto _ : state) {
        auto tets = ChDelaunay3D::CreateTetrahedralization(points);
        benchmark::DoNotOptimize(tets);
    }

    state.SetItemsProcessed(state.iterations() * state.range(0));
}

// Register tests: starting from 50 points up to 5,00
BENCHMARK(BM_Delaunay3D)->RangeMultiplier(10)->Range(50, 5000)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();