"""bioalgs - Bioinformatics algorithms toolkit"""
from typing import Iterator, Iterable

try:
    from importlib.metadata import version
    __version__ = version("bioalgs")
except ImportError:
    # Fallback для старых версий Python
    __version__ = "unknown"

__all__ = [
    "kmergen", "pattern_count", "kmer_count_d", "most_frequent_words",
    "revcomp", "ss_idx", "find_clumps", "find_min_skew_poses", 
    "gen_possible_kmers", "hamming_distance", "kmergen_by_pattern"
]

def kmergen(text: str, kmerlen: int) -> Iterator[str]:
    """
    Генерирует k-меры (подстроки длины kmerlen) из строки text.

    Args:
        text (str): Исходная строка ДНК/РНК.
        kmerlen (int): Длина каждого k-мера.

    Yields:
        Iterator[str]: Последовательные k-меры из строки.
    """
    kmernum: int = len(text) - int(kmerlen) + 1
    for i in range(kmernum):
        yield text[i:i + kmerlen]

def pattern_count(text: str, pattern: str) -> int:
    """
    Подсчитывает количество вхождений подстроки pattern в строке text.

    Args:
        text (str): Строка для поиска.
        pattern (str): Подстрока для поиска и подсчета.

    Returns:
        int: Число вхождений pattern в text.
    """
    kmerlen: int = len(pattern)
    count: int = 0
    for kmer in kmergen(text, kmerlen):
        count += kmer == pattern
    return count

def kmer_count_d(text: str, kmerlen: int) -> dict[str, int]:
    """
    Подсчитывает количество всех уникальных k-меров длины kmerlen в строке text
    без использования collections.Counter.

    Args:
        text (str): Строка ДНК/РНК.
        kmerlen (int): Длина k-мера.

    Returns:
        dict[str, int]: Словарь {k-мер: число вхождений}.
    """
    Counter: dict[str, int] = {}
    for kmer in kmergen(text, kmerlen):
        if kmer in Counter:
            Counter[kmer] += 1
        else:
            Counter[kmer] = 1
    return Counter

def most_frequent_words(d: dict[str, int], sep: str = " ") -> str:
    """
    Возвращает строку наиболее часто встречающихся k-меров (через заданный разделитель),
    если таких несколько — объединяет их в одну строку.

    Args:
        d (dict[str, int]): Словарь {k-мер: частота}.
        sep (str): Разделитель строки (по умолчанию пробел)

    Returns:
        str: Строка всех наиболее частых k-меров, объединённых пробелом или другим разделителем.
    """
    max_len: int = max(d.values())
    words: str = ""
    for k, v in d.items():
        assert v <= max_len
        if v == max_len:
            words += f"{sep}{k}"
    return words.strip(sep)

def revcomp(text: str) -> str:
    """
    Возвращает обратную комплементарную цепь для ДНК-строки text.

    Args:
        text (str): Исходная ДНК-цепь в формате строки.

    Returns:
        str: Обратная комплементарная цепь (реверсированная и заменённая по правилам комплементарности).
    """
    complement_dna = str.maketrans("ATCG", "TAGC")
    return text[::-1].translate(complement_dna)

def ss_idx(ss: str, s: str) -> list[int]:
    """
    Находит все индексы (позиции начала) подстроки ss в строке s.

    Args:
        ss (str): Искомая подстрока.
        s (str): Строка для поиска.

    Returns:
        list[int]: Список индексов вхождений ss в s.
    """
    idx: int = -1
    idx_list: list[int] = []
    while True:
        idx = s.find(ss, idx+1)
        if idx == -1:
            break
        idx_list.append(idx)
    return idx_list

def find_clumps(genome: str, k: int, L: int, t: int) -> list[str]:
    """
    Находит все k-меры, которые встречаются не менее t раз в каждом окне длины L
    в последовательности genome (клампы).

    Args:
        genome (str): Последовательность ДНК.
        k (int): Длина k-мера.
        L (int): Размер окна.
        t (int): Минимальное число повторов.

    Returns:
        list[str]: Список k-меров, обнаруженных в клампах.
    """
    clumps: set[str] = set()  # avoid duplicates
    for i in range(len(genome) - L + 1):
        genomeL = genome[i:i + L]
        kmer_d = kmer_count_d(genomeL, k)
        for kmer, count in kmer_d.items():
            if count >= t:
                clumps.add(kmer)
    return list(clumps)

def find_min_skew_poses(genome: str) -> list[int]:
    """
    Находит позиции минимального skew (разница между количеством G и C по ходу строки),
    обычно применяется для локализации начала репликации в ДНК бактерий.

    Args:
        genome (str): Последовательность ДНК.

    Returns:
        list[int]: Индексы (позиции), где skew минимален.
    """
    min_skew, skew = 0, 0
    poses: list[int] = []
    for k, v in enumerate(genome):
        match v:
            case "G":
                skew += 1
            case "C":
                skew -= 1
                if skew < min_skew:
                    min_skew = skew
                    poses = [k+1]
                elif skew == min_skew:
                    poses.append(k+1)
            case _:
                if skew == min_skew:
                    poses.append(k+1)
    return poses

def gen_possible_kmers(kmerlen: int, bases: str = "ATCG") -> Iterable[str]:
    """
    Генерирует все возможные k-меры заданной длины из алфавита bases.

    Args:
        kmerlen (int): Длина k-мера.
        bases (str, optional): Строка с допустимыми буквами (по умолчанию "ATCG").

    Returns:
        Iterable[str]: Генератор всех возможных k-меров.
    """
    # Есть, возможно, более эффективный вариант с itertools
    # import itertools
    # def gen_possible_kmers(kmerlen: int, bases: str = "ATCG") -> Iterable[str]:
    #     return (''.join(kmer) for kmer in itertools.product(bases, repeat=kmerlen))
    kmers: Iterable[str] = list(bases)
    for _ in range(kmerlen - 1):
        kmers = (kmer + base
                for kmer in kmers
                for base in bases)
    return kmers

def hamming_distance(s1: str, s2: str) -> int:
    """
    Вычисляет расстояние Хэмминга между двумя строками s1 и s2 одинаковой длины,
    используя векторные операции NumPy.

    Args:
        s1 (str): Первая строка.
        s2 (str): Вторая строка.

    Returns:
        int: Число позиций, в которых символы различаются.

    Requires:
        numpy (np) должен быть импортирован.
    """
    import numpy as np # pyright: ignore[reportMissingImports]
    return np.sum(np.frombuffer(s1.encode(), dtype="S1") != np.frombuffer(s2.encode(), dtype="S1")) # type: ignore

def kmergen_by_pattern(_pattern: str = "CACAGC", d: int = 2, bases: str = "ATCG"):
    """
    Генерирует все возможные k-меры с расстоянием Хэмминга не более d от исходного паттерна.
    
    Args:
        _pattern (str): Исходный паттерн для генерации соседних k-меров. По умолчанию "CACAGC".
        d (int): Максимальное расстояние Хэмминга (количество различающихся позиций). По умолчанию 2.
        bases (str): Алфавит допустимых символов. По умолчанию "ATCG".
    
    Returns:
        set[str]: Множество всех k-меров с расстоянием Хэмминга не более d от исходного паттерна.
    """
    _set: set[str] = set()
    _set.add(_pattern)
    for _ in range(d):
        subset: set[str] = set()
        for _pattern in _set:
            for i, char in enumerate(_pattern):
                # print(_pattern[0:i] + _pattern[i].lower() + _pattern[i:])
                other_bases = list(bases)
                other_bases.remove(char)
                # left_set = kmergen_by_pattern(_pattern[0:i])
                # right_set = kmergen_by_pattern(_pattern[i:])
                for base in other_bases:
                    # print(_pattern[0:i] + base.lower() + _pattern[i:])
                    # _set.add(_pattern[0:i] + base + _pattern[i:])
                    #
                    # Думаю, вот здесь надо сделать перебор для оставшихся
                    subset.add(_pattern[0:i] + base + _pattern[i+1:])
        _set.update(subset)
    return _set