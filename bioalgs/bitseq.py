from __future__ import annotations
from typing import Union, Iterator, overload, List, Set
import re

class DNAString:
    code: dict[str, int] = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    inv: dict[int, str]  = {v: k for k, v in code.items()}
    alphalen: int        = len(code)

    def __init__(self, init_value: object) -> None:
        if isinstance(init_value, str):
            self._validate_dna_string(init_value)
            self.text: str = init_value.upper()
            self.num:  int = self._encode(self.text)
        elif isinstance(init_value, int):
            if init_value < 0:
                raise ValueError("Integer must be non-negative")
            self.num = init_value
            self.text = self._decode(self.num)
        else:
            raise TypeError("Argument must be str or int")
        self.len: int = len(self.text)

    @staticmethod
    def _validate_dna_string(s: str) -> None:
        if not re.fullmatch(r'[ACGT]+', s.upper()):
            raise ValueError("Only A,C,G,T allowed")

    def _encode(self, text: str) -> int:
        x = 0
        for ch in text:
            x = (x << 2) | self.code[ch]
        return x

    def _decode(self, num: int) -> str:
        if num == 0:
            return 'A'
        s: list[str] = []
        x = num
        while x:
            s.append(self.inv[x & 3])
            x >>= 2
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

    def __add__(self, other: Union[str, DNAString]) -> DNAString:
        if isinstance(other, DNAString):
            s2 = other.text
        else:
            s2 = other.upper()
        return DNAString(self.text + s2)

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
        
        Алгоритм:
        1. XOR двух закодированных чисел выделяет различающиеся биты
        2. Подсчёт единиц в результате XOR даёт количество различающихся нуклеотидов
        """
        if self.len != other.len:
            raise ValueError("Lengths must match")
        
        # XOR выделяет все различающиеся биты
        # xor_result = self.num ^ other.num
        
        # Подсчитываем количество единиц в результате XOR
        # return bin(xor_result).count('1') // 2
        return (self.num ^ other.num).bit_count() // 2


    def find_motif(self, motif: Union[str, DNAString]) -> List[int]:
        """
        Находит все позиции мотива используя битовые операции и скользящее окно.
        
        Алгоритм:
        1. Кодируем мотив в число
        2. Создаём маску для извлечения подпоследовательности нужной длины
        3. Сдвигаем окно по основной последовательности, сравнивая с мотивом
        """
        if isinstance(motif, DNAString):
            motif_text = motif.text
            motif_num = motif.num
            motif_len = motif.len
        else:
            motif_text = motif.upper()
            motif_len = len(motif_text)
            motif_num = self._encode(motif_text)
        
        if motif_len > self.len:
            return []
        
        positions: List[int] = []
        
        # Создаём маску для извлечения подпоследовательности длины motif_len
        mask = (1 << (2 * motif_len)) - 1 # 1 x 2*motif_len 
        
        # Извлекаем первое окно (начиная с позиции 0)
        shift_amount = 2 * (self.len - motif_len)
        current_window = (self.num >> shift_amount) & mask
        
        # Проверяем первую позицию
        if current_window == motif_num:
            positions.append(0)
        
        # Сдвигаем окно по всей последовательности
        for pos in range(1, self.len - motif_len + 1):
            # Сдвигаем окно на одну позицию вправо
            shift_amount -= 2
            current_window = (self.num >> shift_amount) & mask
            
            if current_window == motif_num:
                positions.append(pos)
        
        return positions


    def neighbors(self, d: int) -> Set[DNAString]:
        encoded = self._get_encoded_neighbors(d)
        return {DNAString(n) for n in encoded}

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
