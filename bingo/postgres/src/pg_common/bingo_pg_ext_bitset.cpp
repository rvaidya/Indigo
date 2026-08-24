#include "bingo_pg_fix_pre.h"

#include "bingo_pg_fix_post.h"

#include "base_c/bitarray.h"
#include "base_cpp/ptr_array.h"
#include "bingo_pg_ext_bitset.h"
#include <algorithm>
#include <climits>

BingoPgExternalBitset::BingoPgExternalBitset()
{
    _initWords(BITS_PER_WORD);
}

BingoPgExternalBitset::BingoPgExternalBitset(int nbits)
{
    _initWords(nbits);
}

void* BingoPgExternalBitset::serialize(int& size)
{
    size = (_length + 2) * sizeof(qword);
    _serializeWords.resize(_length + 2);
    _serializeWords[0] = _bitsNumber;
    _serializeWords[1] = (*_lastWordPtr);
    memcpy(&_serializeWords[2], _words, _length * sizeof(qword));
    return _serializeWords.ptr();
}

void BingoPgExternalBitset::deserialize(void* data_ptr, int data_len, bool ext)
{
    if (data_ptr == nullptr || data_len < (int)(2 * sizeof(qword)))
        throw indigo::Exception("invalid bingo bitset data length %d", data_len);

    qword* data = (qword*)data_ptr;
    if (data[0] == 0 || data[0] > INT_MAX)
        throw indigo::Exception("invalid bingo bitset size %llu", (unsigned long long)data[0]);

    const int bits_number = (int)data[0];
    const int length = _wordIndex(bits_number - 1) + 1;
    const size_t required_size = (size_t)(length + 2) * sizeof(qword);
    if (length <= 0 || required_size > (size_t)data_len)
        throw indigo::Exception("bingo bitset size %d exceeds serialized data length %d", bits_number, data_len);
    if (data[1] > (qword)length)
        throw indigo::Exception("invalid bingo bitset words-in-use value %llu for %d words", (unsigned long long)data[1], length);

    _bitsNumber = bits_number;
    _lastWordPtr = &(data[1]);
    _length = length;
    qword* words_ptr = &(data[2]);
    if (ext)
    {
        _words = words_ptr;
    }
    else
    {
        _internalWords.resize(_length);
        _internalWords.copy(words_ptr, _length);
        _words = _internalWords.ptr();
        _lastWord = data[1];
        _lastWordPtr = &_lastWord;
    }
}

void BingoPgExternalBitset::_recalculateWordsInUse()
{
    // Traverse the bitset until a used word is found
    int i = _length - 1;
    for (; i >= 0; --i)
        if (_words[i] != 0)
            break;
    (*_lastWordPtr) = i + 1; // The new logical size
}

void BingoPgExternalBitset::_initWords(int nbits)
{
    _lastWordPtr = &_lastWord;
    (*_lastWordPtr) = 0;
    _length = _wordIndex(nbits - 1) + 1;
    _internalWords.clear_resize(_length);
    _internalWords.zerofill();
    _bitsNumber = nbits;
    _words = _internalWords.ptr();
}

void BingoPgExternalBitset::_expandTo(int wordIndex)
{
    int wordsRequired = wordIndex + 1;
    if ((*_lastWordPtr) < wordsRequired)
    {
        (*_lastWordPtr) = wordsRequired;
    }
}

int BingoPgExternalBitset::_bitCount(qword b) const
{
    b = (b & 0x5555555555555555LL) + ((b >> 1) & 0x5555555555555555LL);
    b = (b & 0x3333333333333333LL) + ((b >> 2) & 0x3333333333333333LL);
    b = (b + (b >> 4)) & 0x0F0F0F0F0F0F0F0FLL;
    b = b + (b >> 8);
    b = b + (b >> 16);
    b = (b + (b >> 32)) & 0x0000007F;
    return (int)b;
}
