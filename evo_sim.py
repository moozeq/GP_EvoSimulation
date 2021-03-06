#!/usr/bin/env python3
import argparse
import math
import random
import sys
from collections import Counter
from typing import List, Generator

from seq_utils import (transitions, transversions, plot, is_transition, is_transversion, download_sequence,
                       read_sequences, )


def abort(msg: str):
    print(f'[ERROR] {msg}')
    sys.exit(1)


class EvoSimulation:
    DEFAULT_ALPHA = 1.0
    DEFAULT_BETA = 0.25
    DEFAULT_TIME = 1
    DEFAULT_STEPS = 100

    def __init__(self, seq: str, time: float, steps: int):
        self.seq = seq
        self.time = time
        self.steps = steps

    def _sim_JC(self, alpha: float) -> Generator[str, None, None]:
        sim_seq = list(self.seq)
        prob_trans = EvoSimulation._pJC_trans(self.time, alpha)
        for _ in range(self.steps):
            for i, letter in enumerate(sim_seq):
                if random.random() < prob_trans:
                    sub_nucleotides = transitions[letter] + transversions[letter]
                    sub_letter = random.choice(sub_nucleotides)
                    sim_seq[i] = sub_letter
            yield ''.join(sim_seq)

    def _sim_K(self, alpha: float, beta: float) -> Generator[str, None, None]:
        sim_seq = list(self.seq)
        prob_transition = EvoSimulation._pK_transition(self.time, alpha, beta)
        prob_transversion = EvoSimulation._pK_transversion(self.time, alpha, beta)
        for _ in range(self.steps):
            for i, letter in enumerate(sim_seq):
                if (r := random.random()) < prob_transition + prob_transversion:
                    if r < prob_transition:
                        sub_nucleotides = transitions[letter]
                    elif prob_transition <= r < prob_transition + prob_transversion:
                        sub_nucleotides = transversions[letter]
                    else:
                        raise Exception(f'Wrong probability! r = {r}')
                    sub_letter = random.choice(sub_nucleotides)
                    sim_seq[i] = sub_letter
            yield ''.join(sim_seq)

    def run_JC(self, alpha: float, /, *,
               plotting: bool = True, output_file: str = '') -> List[int]:
        print(f'[*] Simulating Jukes-Cantor, alpha = {alpha}')
        dists = [
            EvoSimulation.count_distance(seq)
            for seq in self._sim_JC(alpha)
        ]
        if plotting:
            plot([dists], ['JC69'],
                 title=f'Jukes-Cantor simulation, alpha = {alpha}',
                 xlabel=f'Simulation step, dt = {self.time}',
                 ylabel=f'L1 distance',
                 output_file=output_file)
        return dists

    def run_K(self, alpha: float, beta: float, /, *,
              plotting: bool = True, output_file: str = '') -> List[int]:
        print(f'[*] Simulating Kimura, alpha = {alpha}, beta = {beta}')
        dists = [
            EvoSimulation.count_distance(seq)
            for seq in self._sim_K(alpha, beta)
        ]
        if plotting:
            plot([dists], ['Kimura'],
                 title=f'Kimura simulation, alpha = {alpha}, beta = {beta}',
                 xlabel=f'Simulation step, dt = {self.time}',
                 ylabel=f'L1 distance',
                 output_file=output_file)
        return dists

    @staticmethod
    def count_distance(seq: str) -> int:
        uni_freq = [0.25, 0.25, 0.25, 0.25]
        counter = Counter(seq)
        freqs = (
            counter['A'] / sum(counter.values()),
            counter['C'] / sum(counter.values()),
            counter['T'] / sum(counter.values()),
            counter['G'] / sum(counter.values()),
        )
        l1_dists = sum(
            abs(a - b)
            for a, b in zip(uni_freq, freqs)
        )
        return l1_dists

    @staticmethod
    def _pJC_same(t: float, alpha: float):
        return 0.25 + 0.75 * math.e ** (-4 * alpha * t)

    @staticmethod
    def _pJC_trans(t: float, alpha: float):
        return 0.25 - 0.25 * math.e ** (-4 * alpha * t)

    @staticmethod
    def _pK_same(t: float, alpha: float, beta: float):
        return 0.25 + 0.25 * math.e ** (-4 * beta * t) + 0.5 * math.e ** (-2 * (alpha + beta) * t)

    @staticmethod
    def _pK_transition(t: float, alpha: float, beta: float):
        return 0.25 + 0.25 * math.e ** (-4 * beta * t) - 0.5 * math.e ** (-2 * (alpha + beta) * t)

    @staticmethod
    def _pK_transversion(t: float, alpha: float, beta: float):
        return 0.25 - 0.25 * math.e ** (-4 * beta * t)

    @staticmethod
    def pJC(a: str, b: str, t: float, /,
            alpha: float = DEFAULT_ALPHA):
        same = EvoSimulation._pJC_same(t, alpha)
        trans = EvoSimulation._pJC_trans(t, alpha)

        same_n = 0
        trans_n = 0
        for a1, a2 in zip(a, b):
            if a1 == '-' or a2 == '-':
                continue
            elif a1 == a2:
                same_n += 1
            elif a1 != a2:
                trans_n += 1
            else:
                raise Exception(f'Should not be here, a1 = {a1}, a2 = {a2}')

        probability = (same ** same_n) * (trans ** trans_n)
        return probability

    @staticmethod
    def pK(a: str, b: str, t: float, /,
           alpha: float = DEFAULT_ALPHA, beta: float = DEFAULT_BETA):
        same = EvoSimulation._pK_same(t, alpha, beta)
        transition = EvoSimulation._pK_transition(t, alpha, beta)
        transversion = EvoSimulation._pK_transversion(t, alpha, beta)

        same_n = 0
        transi_n = 0
        transv_n = 0
        for a1, a2 in zip(a, b):
            if a1 == '-' or a2 == '-':
                continue
            elif a1 == a2:
                same_n += 1
            elif is_transition(a1, a2):
                transi_n += 1
            elif is_transversion(a1, a2):
                transv_n += 1
            else:
                raise Exception(f'Should not be here, a1 = {a1}, a2 = {a2}')

        probability = (same ** same_n) * (transition ** transi_n) * (transversion ** transv_n)
        return probability


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evolution models simulator')
    parser.add_argument('seq', type=str, help='sequence fasta filename or NCBI ID')
    parser.add_argument('model', type=str, choices=['jukes-cantor', 'kimura'], help='evolutionary model')
    parser.add_argument('-a', '--alpha', type=float, default=EvoSimulation.DEFAULT_ALPHA, help='alpha param')
    parser.add_argument('-b', '--beta', type=float, default=EvoSimulation.DEFAULT_BETA,
                        help='beta param (for Kimura model)')
    parser.add_argument('-d', '--time', type=float, default=EvoSimulation.DEFAULT_TIME,
                        help='time delta of single step')
    parser.add_argument('-K', '--steps', type=int, default=EvoSimulation.DEFAULT_STEPS,
                        help='how many simulation steps')
    parser.add_argument('-o', '--output', type=str, help='output file')
    args = parser.parse_args()

    if not args.seq.endswith('.fasta'):
        args.seq = download_sequence(args.seq)
    sequences = read_sequences(args.seq)
    output_file = args.output if args.output else ''
    try:
        sim = EvoSimulation(sequences[0], args.time, args.steps)
        if args.model == 'jukes-cantor':
            dists = sim.run_JC(args.alpha, output_file=output_file)
        elif args.model == 'kimura':
            dists = sim.run_K(args.alpha, args.beta, output_file=output_file)
        else:
            raise Exception('Wrong model')
    except Exception as err:
        abort(str(err))
