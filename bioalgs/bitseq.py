from __future__ import annotations
from typing import Union, Iterator, overload, List, Set, Optional
# import re

class DNAString:
    code: dict[str, int] = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    inv: dict[int, str]  = {v: k for k, v in code.items()}
    _valid_bases: set[str] = set(code.keys())
    alphalen: int        = len(code)

    def __init__(self, init_value: object, length: Optional[int] = None) -> None:
        if isinstance(init_value, str):
            if not self.is_valid_sequence(init_value):
                raise ValueError(f"Only {','.join(self.code.keys())} allowed")
            self.text: str = init_value.upper()
            self.len: int = len(self.text)
            self.num: int = self._encode(self.text)
            
            # Проверяем соответствие переданной длины
            if length is not None and length != self.len:
                raise ValueError(f"String length {self.len} doesn't match specified length {length}")
                
        elif isinstance(init_value, int):
            if init_value < 0:
                raise ValueError("Integer must be non-negative")
            if length is None:
                raise ValueError("Length must be specified when creating from integer")
            
            self.num = init_value
            self.len = length
            self.text = self._decode(self.num, self.len)
            
        else:
            raise TypeError("Argument must be str or int")

    @classmethod
    def is_valid_sequence(cls, sequence: str) -> bool:
        """Проверяет, является ли последовательность валидной ДНК"""
        return cls._valid_bases.issuperset(sequence.upper())

    def _encode(self, text: str) -> int:
        x = 0
        for ch in text:
            x = (x << 2) | self.code[ch]
        return x

    def _decode(self, num: int, k: int) -> str:
        """
        Декодирует число в ДНК-строку заданной длины k.
        
        Args:
            num: Закодированное число
            k: Длина результирующей строки
        
        Returns:
            ДНК-строка длины k
        """
        if k == 0:
            return ''
        
        s: list[str] = []
        
        # Извлекаем ровно k символов (справа налево) и переворачиваем
        for _ in range(k):
            ## Забираем два крайних бита и смотрим, что это за буква в словаре
            s.append(self.inv[num & 3]) # 3₁₀ == 11₂ 
            ## Сдвигаемся вправо, убирая прочитанные биты
            num >>= 2                   # 58₁₀ = 00111010₂ >> 2 = 00001110₂ = 14₁₀
        return ''.join(reversed(s))
    
    def __str__(self) -> str:
        return self.text

    def __repr__(self) -> str:
        # return f"DNAString('{self.text}')"
        return self.text

    def __len__(self) -> int:
        return self.len

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return self.text == other.upper()
        if isinstance(other, int):
            return self.num == other
        return False

    def __hash__(self) -> int:
        return hash(self.num)

    @overload
    def __getitem__(self, key: int) -> str: ...
    @overload
    def __getitem__(self, key: slice) -> DNAString: ...
    def __getitem__(self, key: Union[int, slice]) -> Union[str, DNAString]:
        if isinstance(key, int):
            return self.text[key]
        return DNAString(self.text[key])

    def __iter__(self) -> Iterator[str]:
        return iter(self.text)

    def __add__(self, other: Union[int, str, DNAString]) -> DNAString:
        if isinstance(other, int):
            new_num = self.num + other
            new_len = max(self.len, self._min_length_for_number(new_num))
            return DNAString(new_num, new_len)
        
        if isinstance(other, DNAString):
            s2 = other.text
        else:
            s2 = other.upper()
            
        return DNAString(self.text + s2)

    @staticmethod
    def _min_length_for_number(num: int) -> int:
        """Возвращает минимальную длину ДНК строки для представления числа"""
        if num == 0:
            return 1
        
        import math
        return math.ceil(num.bit_length() / 2) # 2 бита на нуклеотид

    def __mul__(self, n: int) -> DNAString:
        if n < 0:
            raise ValueError("Multiplier must be ≥ 0")
        return DNAString(self.text * n)

    def __contains__(self, item: object) -> bool:
        if isinstance(item, DNAString):
            return item.text in self.text
        if isinstance(item, str):
            return item.upper() in self.text
        return False

    # def __lt__(self, other: DNAString) -> bool:
    #     return self.text < other.text

    # def __le__(self, other: DNAString) -> bool:
    #     return self.text <= other.text

    # def __gt__(self, other: DNAString) -> bool:
    #     return self.text > other.text

    # def __ge__(self, other: DNAString) -> bool:
    #     return self.text >= other.text

    # def __bool__(self) -> bool:
    #     return bool(self.text)

    def complement(self) -> DNAString:
        """
        Возвращает комплементарную последовательность. 
        Комплементарность достигается инверсией битов:
        A (00) XOR 11 = T (11)
        C (01) XOR 11 = G (10) 
        G (10) XOR 11 = C (01)
        T (11) XOR 11 = A (00)
        """
        # 1 << n = 2**n
        all_ones = (1 << (2 * self.len)) - 1
        return DNAString(self.num ^ all_ones)


    def reverse_complement(self) -> DNAString:
        return self.complement()[::-1]

    def gc_content(self) -> float:
        return (self.text.count('G') + self.text.count('C')) / self.len if self.len else 0.0

    def hamming_distance(self, other: DNAString) -> int:
        """
        Вычисляет расстояние Хэмминга через битовые операции.
        """
        if self.len != other.len:
            raise ValueError("Lengths must match")
        
        # XOR выделяет все различающиеся биты
        # xor_result = self.num ^ other.num
        
        # Подсчитываем количество единиц в результате XOR
        # return bin(xor_result).count('1') // 2
        
        xor_result = self.num ^ other.num
        distance = 0
        
        # Проверяем каждые 2 бита отдельно
        for _ in range(self.len):
            # Если любой из 2 битов нуклеотида отличается
            if xor_result & 3:  # 3 = 11₂
                distance += 1
            xor_result >>= 2
        
        return distance
    
    def kmergen_nums(self, k: int) -> Iterator[int]:
        """
        Генерирует закодированные k-меры как числа (быстрее для алгоритмов).
        """
        if k <= 0 or k > self.len:
            return
        
        mask = (1 << (2 * k)) - 1
        shift_amount = 2 * (self.len - k)
        
        for _ in range(self.len - k + 1):
            yield (self.num >> shift_amount) & mask
            shift_amount -= 2

    def kmergen(self, k: int) -> Iterator[DNAString]:
        """
        Основной метод генерации k-меров, использующий kmergen_nums(...) для эффективности.
        """
        for kmer_num in self.kmergen_nums(k):
            yield DNAString(kmer_num, k)
    #
    # Ниже разные версии kmergen -- можно удалить после обкатки нового варианта
    # или откатиться
    #
    # def kmergen(self, k: int) -> Set[DNAString]:
    #     return set(
    #         map(DNAString, {self.text[i:i+k] for i in range(len(self.text)-k+1)}))

    # def kmergen(self, k: int) -> Iterator[DNAString]:
    #     """
    #     Генерирует все k-меры данной ДНК последовательности используя битовые операции.
        
    #     Args:
    #     k: Длина k-мера
        
    #     Yields:
    #         DNAString: k-мер в виде объекта DNAString
            
    #     Raises:
    #         ValueError: Если k больше длины последовательности или k <= 0
    #     """
    #     # Валидация входных данных
    #     if k <= 0:
    #         raise ValueError("k must be positive")
    #     if k > self.len:
    #         raise ValueError(f"k ({k}) cannot be greater than sequence length ({self.len})")
        
    #     # Создаём маску для извлечения подпоследовательности длины k
    #     mask = (1 << (2 * k)) - 1
        
    #     # Генерируем все k-меры
    #     for pos in range(self.len - k + 1):
    #         shift_amount = 2 * (self.len - k - pos)
    #         current_window = (self.num >> shift_amount) & mask
    #         yield DNAString(current_window, k)  # Передаём длину!

    #     # Создаём маску для извлечения подпоследовательности длины motif_len
    #     mask = (1 << (2 * k)) - 1 # 1 x 2*k 
        
    #     # Извлекаем первое окно (начиная с позиции 0)
    #     shift_amount = 2 * (self.len - k)
    #     current_window = (self.num >> shift_amount) & mask
        
    #     yield DNAString(current_window)
        
    #     # Сдвигаем окно по всей последовательности
    #     for _ in range(1, self.len - k + 1):
    #         # Сдвигаем окно на одну позицию вправо
    #         shift_amount -= 2
    #         current_window = (self.num >> shift_amount) & mask
    #         yield DNAString(current_window)

    def find_motif(self, motif: Union[str, DNAString]) -> List[int]:
        """
        Находит все позиции мотива в последовательности.
        Использует kmergen для генерации окон.
        """
        if isinstance(motif, DNAString):
            motif_num = motif.num
            motif_len = motif.len
        else:
            motif_text = motif.upper()
            motif_len = len(motif_text)
            motif_num = self._encode(motif_text)
        
        if motif_len > self.len:
            return []
        
        positions: List[int] = []
        
        # Используем kmergen_nums (если реализован) для скорости
        if hasattr(self, 'kmergen_nums'):
            for pos, kmer_num in enumerate(self.kmergen_nums(motif_len)):
                if kmer_num == motif_num:
                    positions.append(pos)
        else:
            # Fallback на обычный kmergen
            for pos, kmer in enumerate(self.kmergen(motif_len)):
                if kmer.num == motif_num:
                    positions.append(pos)
        
        return positions

    def neighbors(self, d: int) -> Set[DNAString]:
        encoded = self._get_encoded_neighbors(d)
        return {DNAString(n, self.len) for n in encoded}

    def _get_encoded_neighbors(self, d: int) -> Set[int]:
        current: Set[int] = {self.num}
        for _ in range(d):
            next_set: Set[int] = set()
            for val in current:
                for pos in range(self.len):
                    shift = 2 * (self.len - 1 - pos)
                    orig = (val >> shift) & 3
                    for b in range(self.alphalen):
                        if b == orig:
                            continue
                        mask = ~(3 << shift)
                        new = (val & mask) | (b << shift)
                        next_set.add(new)
            current |= next_set
        return current
    
    
