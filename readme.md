
# Exparser

A Python library for the analysis of experimental data.

Copyright 2011-2014 Sebastiaan Mathôt

## Overview


- [class __DataMatrix__](#class-__datamatrix__)
	- [function __DataMatrix\.\_\_add\_\___\(dm, cautious=False\)](#function-__datamatrix__add____dm-cautiousfalse)
	- [function __DataMatrix\.\_\_getitem\_\___\(key\)](#function-__datamatrix__getitem____key)
	- [function __DataMatrix\.\_\_init\_\___\(a, structured=False\)](#function-__datamatrix__init____a-structuredfalse)
	- [function __DataMatrix\.\_\_iter\_\___\(\)](#function-__datamatrix__iter____)
	- [function __DataMatrix\.\_\_len\_\___\(\)](#function-__datamatrix__len____)
	- [function __DataMatrix\.\_\_setitem\_\___\(key, val\)](#function-__datamatrix__setitem____key-val)
	- [function __DataMatrix\.addField__\(key, default=None, dtype=<type 'numpy\.int32'>\)](#function-__datamatrixaddfield__key-defaultnone-dtypetype-numpyint32)
	- [function __DataMatrix\.asArray__\(\)](#function-__datamatrixasarray__)
	- [function __DataMatrix\.balance__\(key, maxErr, ref=0, verbose=False\)](#function-__datamatrixbalance__key-maxerr-ref0-verbosefalse)
	- [function __DataMatrix\.calcPerc__\(vName, targetVName, keys=None, nBin=None\)](#function-__datamatrixcalcperc__vname-targetvname-keysnone-nbinnone)
	- [function __DataMatrix\.collapse__\(keys, vName\)](#function-__datamatrixcollapse__keys-vname)
	- [function __DataMatrix\.columns__\(dtype=False\)](#function-__datamatrixcolumns__dtypefalse)
	- [function __DataMatrix\.count__\(key\)](#function-__datamatrixcount__key)
	- [function __DataMatrix\.explode__\(N\)](#function-__datamatrixexplode__n)
	- [function __DataMatrix\.group__\(keys, \_sort=True\)](#function-__datamatrixgroup__keys-_sorttrue)
	- [function __DataMatrix\.intertrialer__\(keys, dv, \_range=\[1\]\)](#function-__datamatrixintertrialer__keys-dv-_range1)
	- [function __DataMatrix\.range__\(\)](#function-__datamatrixrange__)
	- [function __DataMatrix\.recode__\(key, coding\)](#function-__datamatrixrecode__key-coding)
	- [function __DataMatrix\.removeField__\(key\)](#function-__datamatrixremovefield__key)
	- [function __DataMatrix\.removeNan__\(key\)](#function-__datamatrixremovenan__key)
	- [function __DataMatrix\.rename__\(oldKey, newKey\)](#function-__datamatrixrename__oldkey-newkey)
	- [function __DataMatrix\.select__\(query, verbose=True\)](#function-__datamatrixselect__query-verbosetrue)
	- [function __DataMatrix\.selectByStdDev__\(keys, dv, thr=2\.5, verbose=False\)](#function-__datamatrixselectbystddev__keys-dv-thr25-verbosefalse)
	- [function __DataMatrix\.selectColumns__\(keys\)](#function-__datamatrixselectcolumns__keys)
	- [function __DataMatrix\.shuffle__\(\)](#function-__datamatrixshuffle__)
	- [function __DataMatrix\.sort__\(keys, ascending=True\)](#function-__datamatrixsort__keys-ascendingtrue)
	- [function __DataMatrix\.split__\(key\)](#function-__datamatrixsplit__key)
	- [function __DataMatrix\.ttest__\(keys, dv, paired=True, collapse=None\)](#function-__datamatrixttest__keys-dv-pairedtrue-collapsenone)
	- [function __DataMatrix\.unique__\(key\)](#function-__datamatrixunique__key)
	- [function __DataMatrix\.where__\(query\)](#function-__datamatrixwhere__query)
	- [function __DataMatrix\.withinize__\(vName, targetVName, key, verbose=True, whiten=False\)](#function-__datamatrixwithinize__vname-targetvname-key-verbosetrue-whitenfalse)



<span class="ClassDoc YAMLDoc" id="DataMatrix" markdown="1">

# class __DataMatrix__

Provides functionality for convenient processing of experimental data.

<span class="FunctionDoc YAMLDoc" id="DataMatrix-__add__" markdown="1">

## function __DataMatrix\.\_\_add\_\___\(dm, cautious=False\)

Concatenates two DataMatrices. Implements the + operator.

__Example:__

~~~ .python
dm3 = dm1 + dm2
~~~

__Arguments:__

- `dm` -- The DataMatrix to be appended.
	- Type: DataMatrix

__Keywords:__

- `cautious` -- DEPRECATED
	- Default: False

__Returns:__

The concatenation of the current and the passed DataMatrix.

- Type: DataMatrix

</span>

[DataMatrix.__add__]: #DataMatrix-__add__
[__add__]: #DataMatrix-__add__

<span class="FunctionDoc YAMLDoc" id="DataMatrix-__getitem__" markdown="1">

## function __DataMatrix\.\_\_getitem\_\___\(key\)

Returns a column, index, or slice. Note that some operations return a copy of the DataMatrix, so they cannot be used to modify the contents of the DataMatrix.

__Example:__

~~~ .python
dm['rt'] # returns column 'rt' as numpy array
dm[0] # returns first row as DataMatrix
dm[0:2] # returns first two rows as DataMatrix
dm[0]['rt'] = 1 # This doesn't alter the original DataMatrix
dm['rt'][0] = 1 # This does!
~~~

__Arguments:__

- `key` -- A column name or index.
	- Type: int, str, unicode

__Returns:__

A DataMatrix or NumPy array corresponding to a slice or column from this DataMatrix.

- Type: DataMatrix, ndarray

</span>

[DataMatrix.__getitem__]: #DataMatrix-__getitem__
[__getitem__]: #DataMatrix-__getitem__

<span class="FunctionDoc YAMLDoc" id="DataMatrix-__init__" markdown="1">

## function __DataMatrix\.\_\_init\_\___\(a, structured=False\)

Creates a new DataMatrix object.

__Arguments:__

- `a` -- A NumPy array, list, or filename. For unstructured NumPy arrays or lists, the first row is assumed to contain column names. Filenames are assumed to refer to a `.npy` file.
	- Type: ndarray, list, str, unicode

__Keywords:__

- `structured` -- Indicates whether `a` is a structured NumPy array.
	- Default: False
	- Type: bool

</span>

[DataMatrix.__init__]: #DataMatrix-__init__
[__init__]: #DataMatrix-__init__

<span class="FunctionDoc YAMLDoc" id="DataMatrix-__iter__" markdown="1">

## function __DataMatrix\.\_\_iter\_\___\(\)

Implements an iterator for 'for' loops to walk through a DataMatrix row by row.

__Example:__

~~~ .python
for rowDm in dm:
        print rowDm
~~~

__Returns:__

No description

- Type: DataMatrixIterator

</span>

[DataMatrix.__iter__]: #DataMatrix-__iter__
[__iter__]: #DataMatrix-__iter__

<span class="FunctionDoc YAMLDoc" id="DataMatrix-__len__" markdown="1">

## function __DataMatrix\.\_\_len\_\___\(\)

No description specified.

__Returns:__

The number of rows.

- Type: int

</span>

[DataMatrix.__len__]: #DataMatrix-__len__
[__len__]: #DataMatrix-__len__

<span class="FunctionDoc YAMLDoc" id="DataMatrix-__setitem__" markdown="1">

## function __DataMatrix\.\_\_setitem\_\___\(key, val\)

Set a certain variable. Implements assigment.

__Example:__

~~~ .python
dm['rt'] = 100
~~~

__Arguments:__

- `key` -- The name of a key.
	- Type: str, unicode
- `val` -- An array with the new values, or a single new value to use for the entire column.

</span>

[DataMatrix.__setitem__]: #DataMatrix-__setitem__
[__setitem__]: #DataMatrix-__setitem__

<span class="FunctionDoc YAMLDoc" id="DataMatrix-addField" markdown="1">

## function __DataMatrix\.addField__\(key, default=None, dtype=<type 'numpy\.int32'>\)

Creates a new DataMatrix that is a copy of the current DataMatrix with an additional field.

__Example:__

~~~ .python
dm = dm.addField('rt', dtype=float, default=-1000)
~~~

__Source(s):__

- <h>
- <t>
- <t>
- <p>
- <:>
- </>
- </>
- <s>
- <t>
- <a>
- <c>
- <k>
- <o>
- <v>
- <e>
- <r>
- <f>
- <l>
- <o>
- <w>
- <.>
- <c>
- <o>
- <m>
- </>
- <q>
- <u>
- <e>
- <s>
- <t>
- <i>
- <o>
- <n>
- <s>
- </>
- <1>
- <2>
- <0>
- <1>
- <8>
- <1>
- <7>
- </>
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- < >
- <a>
- <d>
- <d>
- <i>
- <n>
- <g>
- <->
- <a>
- <->
- <f>
- <i>
- <e>
- <l>
- <d>
- <->
- <t>
- <o>
- <->
- <a>
- <->
- <s>
- <t>
- <r>
- <u>
- <c>
- <t>
- <u>
- <r>
- <e>
- <d>
- <->
- <n>
- <u>
- <m>
- <p>
- <y>
- <->
- <a>
- <r>
- <r>
- <a>
- <y>
- <>>
__Arguments:__

- `key` -- The name of the new field.
	- Type: str, unicode

__Keywords:__

- `default` -- The default value or `None` for no default.
	- Default: None
- `dtype` -- The dtype for the new field.
	- Default: <type 'numpy.int32'>

__Returns:__

No description

- Type: DataMatrix.

</span>

[DataMatrix.addField]: #DataMatrix-addField
[addField]: #DataMatrix-addField

<span class="FunctionDoc YAMLDoc" id="DataMatrix-asArray" markdown="1">

## function __DataMatrix\.asArray__\(\)

No description specified.

__Returns:__

An array representation of the current DataMatrix.

- Type: ndarray

</span>

[DataMatrix.asArray]: #DataMatrix-asArray
[asArray]: #DataMatrix-asArray

<span class="FunctionDoc YAMLDoc" id="DataMatrix-balance" markdown="1">

## function __DataMatrix\.balance__\(key, maxErr, ref=0, verbose=False\)

Filters the data such that a given column is on average close to a reference value, and is symetrically distributed.

__Arguments:__

- `maxErr` -- The maximum mean error relative to the reference value.
	- Type: int, float
- `key` -- The key to balance. It must have a numeric dtype.
	- Type: str, unicode

__Keywords:__

- `ref` -- The reference value.
	- Default: 0
	- Type: int, float
- `verbose` -- Indicates whether verbose output is printed.
	- Default: False
	- Type: bool

__Returns:__

A balanced copy of the current DataMatrix.

- Type: DataMatrix

</span>

[DataMatrix.balance]: #DataMatrix-balance
[balance]: #DataMatrix-balance

<span class="FunctionDoc YAMLDoc" id="DataMatrix-calcPerc" markdown="1">

## function __DataMatrix\.calcPerc__\(vName, targetVName, keys=None, nBin=None\)

Calculates percentile scores for a variable.

__Arguments:__

- `targetVName` -- The variable to store the percentile scores in. This variable must exist, it is not created.
	- Type: str, unicode
- `vName` -- The variable to calculate percentile scores for.
	- Type: str, unicode

__Keywords:__

- `keys` -- A key or list of keys to split the data by, before calculating percentile scores, so you can calculate scores individually per subject, condition, etc.
	- Default: None
	- Type: list, str, unicode
- `nBin` -- The number of bins or None for continuous scores.
	- Default: None
	- Type: int, NoneType

__Returns:__

No description

- Type: DataMatrix

</span>

[DataMatrix.calcPerc]: #DataMatrix-calcPerc
[calcPerc]: #DataMatrix-calcPerc

<span class="FunctionDoc YAMLDoc" id="DataMatrix-collapse" markdown="1">

## function __DataMatrix\.collapse__\(keys, vName\)

Collapse the data by a (list of) keys and get statistics on a dependent variable.

__Arguments:__

- `keys` -- A key or list of keys to collapse the data on.
	- Type: list, str, unicode
- `vName` -- The dependent variable to collapse. Alternative, you can specifiy a function, in which case the error will be 0.
	- Type: str, unicode, function

__Returns:__

A DataMatrix with the collapsed data, with the descriptives statistics on `vName`.

- Type: DataMatrix

</span>

[DataMatrix.collapse]: #DataMatrix-collapse
[collapse]: #DataMatrix-collapse

<span class="FunctionDoc YAMLDoc" id="DataMatrix-columns" markdown="1">

## function __DataMatrix\.columns__\(dtype=False\)

Returns a list of the columns.

__Keywords:__

- `dtype` -- Indicates whether the dtype for each column should be returned as well.
	- Default: False
	- Type: bool

__Returns:__

If dtype == False: A list of names
If dtype == True: A list of (name, dtype) tuples

- Type: list

</span>

[DataMatrix.columns]: #DataMatrix-columns
[columns]: #DataMatrix-columns

<span class="FunctionDoc YAMLDoc" id="DataMatrix-count" markdown="1">

## function __DataMatrix\.count__\(key\)

Returns the number of different values for a given variable.

__Arguments:__

- `key` -- The variable to count.
	- Type: str, unicode

__Returns:__

The number of different values for 'key'.

- Type: int

</span>

[DataMatrix.count]: #DataMatrix-count
[count]: #DataMatrix-count

<span class="FunctionDoc YAMLDoc" id="DataMatrix-explode" markdown="1">

## function __DataMatrix\.explode__\(N\)

Break up the DataMatrix in N smaller DataMatrices. For splitting a DataMatrix based on column values, see [DataMatrix.split].

__Arguments:__

- `N` -- The number of DataMatrices to explode in.
	- Type: int

__Returns:__

A list of DataMatrices.

- Type: list

</span>

[DataMatrix.explode]: #DataMatrix-explode
[explode]: #DataMatrix-explode

<span class="FunctionDoc YAMLDoc" id="DataMatrix-group" markdown="1">

## function __DataMatrix\.group__\(keys, \_sort=True\)

Split the data into different groups based on unique values for the key variables.

__Arguments:__

- `keys` -- A key or list of keys to split the data on.
	- Type: str, unicode, list

__Keywords:__

- `_sort` -- Indicates whether the groups should be sorted by values.
	- Default: True
	- Type: bool

__Returns:__

A list of DataMatrices.

- Type: list

</span>

[DataMatrix.group]: #DataMatrix-group
[group]: #DataMatrix-group

<span class="FunctionDoc YAMLDoc" id="DataMatrix-intertrialer" markdown="1">

## function __DataMatrix\.intertrialer__\(keys, dv, \_range=\[1\]\)

Adds columns that contain values from the previous or next trial. These columns are called '[dv]_p1' for the next value, '[dv]_m1' for the previous one, etc.

__Arguments:__

- `keys` -- A key or list of keys that define the trial order.
	- Type: list, str, unicode
- `dv` -- The dependent variable.
	- Type: str, unicode

__Keywords:__

- `_range` -- A list of integers that specifies the range for which the operation should be executed.
	- Default: [1]
	- Type: list

__Returns:__

No description

- Type: DataMatrix

</span>

[DataMatrix.intertrialer]: #DataMatrix-intertrialer
[intertrialer]: #DataMatrix-intertrialer

<span class="FunctionDoc YAMLDoc" id="DataMatrix-range" markdown="1">

## function __DataMatrix\.range__\(\)

Gives a list of indices to walk through the current DataMatrix.

__Returns:__

A list of indices.

</span>

[DataMatrix.range]: #DataMatrix-range
[range]: #DataMatrix-range

<span class="FunctionDoc YAMLDoc" id="DataMatrix-recode" markdown="1">

## function __DataMatrix\.recode__\(key, coding\)

Recodes values (i.e. changes one value into another for a given set of columns).

__Arguments:__

- `coding` -- An (oldValue, newValue) tuple, a list of tuples to handle multiple recodings in one go, or a function that takes a value and returns the recoded value.
	- Type: tuple, list, function
- `key` -- The name of the variable to recode, or a list of names to recode multiple variables in one go.
	- Type: str, unicode, list

</span>

[DataMatrix.recode]: #DataMatrix-recode
[recode]: #DataMatrix-recode

<span class="FunctionDoc YAMLDoc" id="DataMatrix-removeField" markdown="1">

## function __DataMatrix\.removeField__\(key\)

Return a DataMatrix that is a copy of the current DataMatrix without the specified field.

__Arguments:__

- `key` -- The name of the field to be removed.
	- Type: str, unicode

__Returns:__

No description

- Type: DataMatrix

</span>

[DataMatrix.removeField]: #DataMatrix-removeField
[removeField]: #DataMatrix-removeField

<span class="FunctionDoc YAMLDoc" id="DataMatrix-removeNan" markdown="1">

## function __DataMatrix\.removeNan__\(key\)

Remove all rows where the specified key is nan.

__Arguments:__

- `key` -- A key that should not have any nan values.
	- Type: str, unicode

__Returns:__

No description

- Type: DataMatrix

</span>

[DataMatrix.removeNan]: #DataMatrix-removeNan
[removeNan]: #DataMatrix-removeNan

<span class="FunctionDoc YAMLDoc" id="DataMatrix-rename" markdown="1">

## function __DataMatrix\.rename__\(oldKey, newKey\)

Renames a column. This function operates in place, so it modifies the current dataMatrix.

__Arguments:__

- `newKey` -- The new name.
	- Type: str, unicode
- `oldKey` -- The old name.
	- Type: str, unicode

__Returns:__

No description

- Type: DataMatrix

</span>

[DataMatrix.rename]: #DataMatrix-rename
[rename]: #DataMatrix-rename

<span class="FunctionDoc YAMLDoc" id="DataMatrix-select" markdown="1">

## function __DataMatrix\.select__\(query, verbose=True\)

Select a subset of the data.

__Arguments:__

- `query` -- A query, e.g. 'rt > 1000'.
	- Type: str, unicode

__Keywords:__

- `verbose` -- Indicates if a summary should be printed.
	- Default: True
	- Type: bool

__Returns:__

No description

- Type: DataMatrix

</span>

[DataMatrix.select]: #DataMatrix-select
[select]: #DataMatrix-select

<span class="FunctionDoc YAMLDoc" id="DataMatrix-selectByStdDev" markdown="1">

## function __DataMatrix\.selectByStdDev__\(keys, dv, thr=2\.5, verbose=False\)

Select only those rows where the value of a given column is within a certain distance from the mean.

__Arguments:__

- `keys` -- A list of keys to create groups for which the deviation is calculated seperately.
	- Type: list
- `dv` -- The dependent variable.
	- Type: str, unicode

__Keywords:__

- `thr` -- The stddev threshold.
	- Default: 2.5
	- Type: float, int
- `verbose` -- Indicates whether detailed output should be provided.
	- Default: False
	- Type: bool

__Returns:__

No description

- Type: DataMatrix

</span>

[DataMatrix.selectByStdDev]: #DataMatrix-selectByStdDev
[selectByStdDev]: #DataMatrix-selectByStdDev

<span class="FunctionDoc YAMLDoc" id="DataMatrix-selectColumns" markdown="1">

## function __DataMatrix\.selectColumns__\(keys\)

Creates a new DataMatrix with only the specified columns.

__Arguments:__

- `keys` -- A column or list of columns to select.
	- Type: list, str, unicode

__Returns:__

No description

- Type: DataMatrix

</span>

[DataMatrix.selectColumns]: #DataMatrix-selectColumns
[selectColumns]: #DataMatrix-selectColumns

<span class="FunctionDoc YAMLDoc" id="DataMatrix-shuffle" markdown="1">

## function __DataMatrix\.shuffle__\(\)

Shuffles the DataMatrix in place.

</span>

[DataMatrix.shuffle]: #DataMatrix-shuffle
[shuffle]: #DataMatrix-shuffle

<span class="FunctionDoc YAMLDoc" id="DataMatrix-sort" markdown="1">

## function __DataMatrix\.sort__\(keys, ascending=True\)

Sorts the DataMatrix in place.

__Arguments:__

- `keys` -- A key or list of keys to use for sorting. The first key is dominant, the second key is next-to-dominant, etc.
	- Type: str, unicode, list

__Keywords:__

- `ascending` -- Indicates whether the sorting should occur in ascending (True) or descending (False) order.
	- Default: True
	- Type: bool

</span>

[DataMatrix.sort]: #DataMatrix-sort
[sort]: #DataMatrix-sort

<span class="FunctionDoc YAMLDoc" id="DataMatrix-split" markdown="1">

## function __DataMatrix\.split__\(key\)

Splits the DataMatrix in chunks such that each chunk only has the same value for the specified column. For splitting a DataMatrix into equally sized parts, see [DataMatrix.explode].

__Arguments:__

- `key` -- A key to split by.
	- Type: str, unicode

__Returns:__

A list of DataMatrices.

- Type: list

</span>

[DataMatrix.split]: #DataMatrix-split
[split]: #DataMatrix-split

<span class="FunctionDoc YAMLDoc" id="DataMatrix-ttest" markdown="1">

## function __DataMatrix\.ttest__\(keys, dv, paired=True, collapse=None\)

Performs t-tests between groups defined by a list of keys.

__Arguments:__

- `keys` -- A list of keys to define the groups.
	- Type: list
- `dv` -- The dependent variable.
	- Type: str, unicode

__Keywords:__

- `paired` -- Determines whether a paired-samples t-test or an independent samples t-test should be conducted.
	- Default: True
	- Type: bool
- `collapse` -- A key to collapse the data on, so that you can do t-tests on (subject) means.
	- Default: None
	- Type: str, unicode, NoneType

__Returns:__

A list of (desc, t, p) tuples.

- Type: list

</span>

[DataMatrix.ttest]: #DataMatrix-ttest
[ttest]: #DataMatrix-ttest

<span class="FunctionDoc YAMLDoc" id="DataMatrix-unique" markdown="1">

## function __DataMatrix\.unique__\(key\)

Gives all unique values for a particular key.

__Arguments:__

- `key` -- A column name.
	- Type: str, unicode

__Returns:__

A list of unique values.

- Type: list

</span>

[DataMatrix.unique]: #DataMatrix-unique
[unique]: #DataMatrix-unique

<span class="FunctionDoc YAMLDoc" id="DataMatrix-where" markdown="1">

## function __DataMatrix\.where__\(query\)

Return indices corresponding to the query.

__Arguments:__

- `query` -- A query, e.g. 'rt > 1000'.
	- Type: str, unicode

__Returns:__

Indices.

- Type: ndarray

</span>

[DataMatrix.where]: #DataMatrix-where
[where]: #DataMatrix-where

<span class="FunctionDoc YAMLDoc" id="DataMatrix-withinize" markdown="1">

## function __DataMatrix\.withinize__\(vName, targetVName, key, verbose=True, whiten=False\)

Removes the between factor variance for a given key (such as subject or file) for a given dependent variable. This operation acts in place.

__Arguments:__

- `targetVName` -- The target variable to store withinized values. This variable should exist.
	- Type: str, unicode
- `vName` -- The dependent variable to withinize.
	- Type: str, unicode
- `key` -- The key that defines the group.
	- Type: str, unicode

__Keywords:__

- `verbose` -- Indicates whether the results should be printed.
	- Default: True
	- Type: bool
- `whiten` -- Indicates whether the data should be whitened so that the standard deviation is 1 and the mean 0.
	- Default: False
	- Type: bool

__Returns:__

No description

- Type: DataMatrix

</span>

[DataMatrix.withinize]: #DataMatrix-withinize
[withinize]: #DataMatrix-withinize

</span>

[DataMatrix]: #DataMatrix




[Exparser]: #exparser
[Overview]: #overview
[class __DataMatrix__]: #class-__datamatrix__
[function __DataMatrix\.\_\_add\_\___\(dm, cautious=False\)]: #function-__datamatrix__add____dm-cautiousfalse
[function __DataMatrix\.\_\_getitem\_\___\(key\)]: #function-__datamatrix__getitem____key
[function __DataMatrix\.\_\_init\_\___\(a, structured=False\)]: #function-__datamatrix__init____a-structuredfalse
[function __DataMatrix\.\_\_iter\_\___\(\)]: #function-__datamatrix__iter____
[function __DataMatrix\.\_\_len\_\___\(\)]: #function-__datamatrix__len____
[function __DataMatrix\.\_\_setitem\_\___\(key, val\)]: #function-__datamatrix__setitem____key-val
[function __DataMatrix\.addField__\(key, default=None, dtype=<type 'numpy\.int32'>\)]: #function-__datamatrixaddfield__key-defaultnone-dtypetype-numpyint32
[function __DataMatrix\.asArray__\(\)]: #function-__datamatrixasarray__
[function __DataMatrix\.balance__\(key, maxErr, ref=0, verbose=False\)]: #function-__datamatrixbalance__key-maxerr-ref0-verbosefalse
[function __DataMatrix\.calcPerc__\(vName, targetVName, keys=None, nBin=None\)]: #function-__datamatrixcalcperc__vname-targetvname-keysnone-nbinnone
[function __DataMatrix\.collapse__\(keys, vName\)]: #function-__datamatrixcollapse__keys-vname
[function __DataMatrix\.columns__\(dtype=False\)]: #function-__datamatrixcolumns__dtypefalse
[function __DataMatrix\.count__\(key\)]: #function-__datamatrixcount__key
[function __DataMatrix\.explode__\(N\)]: #function-__datamatrixexplode__n
[function __DataMatrix\.group__\(keys, \_sort=True\)]: #function-__datamatrixgroup__keys-_sorttrue
[function __DataMatrix\.intertrialer__\(keys, dv, \_range=\[1\]\)]: #function-__datamatrixintertrialer__keys-dv-_range1
[function __DataMatrix\.range__\(\)]: #function-__datamatrixrange__
[function __DataMatrix\.recode__\(key, coding\)]: #function-__datamatrixrecode__key-coding
[function __DataMatrix\.removeField__\(key\)]: #function-__datamatrixremovefield__key
[function __DataMatrix\.removeNan__\(key\)]: #function-__datamatrixremovenan__key
[function __DataMatrix\.rename__\(oldKey, newKey\)]: #function-__datamatrixrename__oldkey-newkey
[function __DataMatrix\.select__\(query, verbose=True\)]: #function-__datamatrixselect__query-verbosetrue
[function __DataMatrix\.selectByStdDev__\(keys, dv, thr=2\.5, verbose=False\)]: #function-__datamatrixselectbystddev__keys-dv-thr25-verbosefalse
[function __DataMatrix\.selectColumns__\(keys\)]: #function-__datamatrixselectcolumns__keys
[function __DataMatrix\.shuffle__\(\)]: #function-__datamatrixshuffle__
[function __DataMatrix\.sort__\(keys, ascending=True\)]: #function-__datamatrixsort__keys-ascendingtrue
[function __DataMatrix\.split__\(key\)]: #function-__datamatrixsplit__key
[function __DataMatrix\.ttest__\(keys, dv, paired=True, collapse=None\)]: #function-__datamatrixttest__keys-dv-pairedtrue-collapsenone
[function __DataMatrix\.unique__\(key\)]: #function-__datamatrixunique__key
[function __DataMatrix\.where__\(query\)]: #function-__datamatrixwhere__query
[function __DataMatrix\.withinize__\(vName, targetVName, key, verbose=True, whiten=False\)]: #function-__datamatrixwithinize__vname-targetvname-key-verbosetrue-whitenfalse