#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######## GENERIC UTILITARIES ############
def applyMappingToColumn(df, column, mapping):
    df[column] = df[column].astype(int)
    df[column] = df.apply( lambda x: mapping[x[column]], axis = 1)
    return df
