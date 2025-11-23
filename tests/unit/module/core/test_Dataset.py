import pytest
import pandas as pd
import numpy as np
from module.core.Dataset import (
    CustomDataFrame,
    DataframeWrapperMixin,
    SelectionError,
    CachedDataFrame,
)
from typing import List, Dict, Any


@pytest.fixture
def sample_data() -> Dict[str, List[Any]]:
    """Create sample data with various data types for testing."""
    return {
        "id": [1, 2, 3, 4, 5],
        "name": ["Alice", "Bob", "Charlie", None, "Eve"],
        "age": [25, 30, 35, 40, None],
        "score": [85.5, 92.0, 78.5, None, 95.0],
        "active": [True, True, False, True, False],
        "category": pd.Categorical(
            ["A", "B", "A", "C", "B"], categories=["A", "B", "C"], ordered=False
        ),
    }


@pytest.fixture
def sample_df(sample_data) -> CustomDataFrame:
    """Create a sample CustomDataFrame for testing."""
    return CustomDataFrame(sample_data)


class TestCustomDataFrame:
    def test_select_literal(self, sample_df):
        result = sample_df.select(active=True)
        assert len(result) == 3
        assert all(result["active"] == True)

    def test_select_iterable(self, sample_df):
        result = sample_df.select(name=["Alice", "Bob"])
        assert len(result) == 2
        assert set(result["name"]) == {"Alice", "Bob"}

    def test_select_na(self, sample_df):
        result = sample_df.select(name="na")
        assert len(result) == 1
        assert pd.isna(result["name"].iloc[0])

    def test_select_notna(self, sample_df):
        result = sample_df.select(score="notna")
        assert len(result) == 4
        assert all(pd.notna(result["score"]))

    def test_select_callable(self, sample_df):
        result = sample_df.select(age=lambda x: x > 30)
        assert len(result) == 2
        assert all(result["age"] > 30)

    def test_select_one(self, sample_df):
        result = sample_df.select(select_one=True, id=1)
        assert isinstance(result, pd.Series)
        assert result["name"] == "Alice"

    def test_select_empty_no_error(self, sample_df):
        result = sample_df.select(select_one=False, id=6)
        assert result.empty

    def test_select_one_error_multiple(self, sample_df):
        with pytest.raises(SelectionError):
            sample_df.select(select_one=True, active=True)

    def test_select_one_error_none(self, sample_df):
        with pytest.raises(SelectionError):
            sample_df.select(select_one=True, name="Nonexistent")

    def test_select_categorical_preservation(self, sample_df):
        result = sample_df.select(category="A")
        assert pd.api.types.is_categorical_dtype(result["category"])
        assert set(result["category"].cat.categories) == {"A"}

    def test_extend_common_columns(self, sample_df):
        other_data = {
            "id": [1, 2, 6],
            "department": ["HR", "Engineering", "Marketing"],
            "salary": [50000, 60000, 55000],
        }
        other_df = CustomDataFrame(other_data)
        result = sample_df.extend(other_df)
        assert "department" in result.columns
        assert result["department"].iloc[0] == "HR"

    def test_extend_no_common_columns(self, sample_df):
        other_data = {"department": ["HR", "Engineering"], "location": ["NY", "SF"]}
        other_df = CustomDataFrame(other_data)
        result = sample_df.extend(other_df)
        assert len(result) == len(sample_df) * len(other_df)
        assert "department" in result.columns
        assert "location" in result.columns


class TestDataFrameWrapperMixin:
    class DataframeWrapperMixinSubclass(DataframeWrapperMixin):
        """Concrete implementation for testing the mixin."""

        def __init__(self, df):
            self._df = df

        @property
        def df(self) -> CustomDataFrame:
            return self._df

    @pytest.fixture
    def wrapper(self, sample_df) -> DataframeWrapperMixinSubclass:
        return self.DataframeWrapperMixinSubclass(sample_df)

    def test_select_through_wrapper(self, wrapper, sample_df):
        result = wrapper.select(active=True)
        expected = sample_df.select(active=True)
        pd.testing.assert_frame_equal(result, expected)

    def test_extend_through_wrapper(self, wrapper):
        other_data = {"id": [1, 2], "department": ["HR", "Engineering"]}
        other_df = CustomDataFrame(other_data)
        result = wrapper.extend(other_df)
        assert len(result) == 2
        assert "department" in result.columns

    def test_list_property(self, wrapper, sample_data):
        result = wrapper.list
        assert isinstance(result, list)
        assert len(result) == len(sample_data["id"])
        assert all(isinstance(item, dict) for item in result)

    def test_contains_operator(self, wrapper):
        assert "name" in wrapper
        assert "nonexistent" not in wrapper


class TestEdgeCases:
    def test_empty_dataframe(self):
        empty_df = CustomDataFrame()
        with pytest.raises(SelectionError):
            empty_df.select(select_one=True)

    def test_all_none_values(self):
        df = CustomDataFrame({"col1": [None, None, None], "col2": [1, 2, 3]})
        result = df.select(col1="na")
        assert len(result) == 3

    def test_extend_with_duplicate_columns(self):
        df1 = CustomDataFrame({"id": [1, 2], "measure": [10, 20]})
        df2 = CustomDataFrame({"id": [1, 3], "region": [100, 300]})
        result = df1.extend(df2)
        assert result.columns.tolist() == ["id", "measure", "region"]
