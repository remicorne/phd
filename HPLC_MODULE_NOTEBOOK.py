#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######## GENERIC UTILITARIES ############
def applyMappingToColumn(df, column, mapping):
    df[column] = int(df[column])
    df[column] = df.apply( lambda x: mapping[x[column]], axis = 1)
    return df
