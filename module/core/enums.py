from enum import StrEnum


class DatasetColumn(StrEnum):
    VALUE = "value"
    UNIT = "unit"
    SUBJECT_ID = "subject_id"
    GROUP_ID = "group_id"
    GROUP_NAME = "group_name"

    @classmethod
    def get_mandatory_columns(cls):
        return [cls.SUBJECT_ID.value, cls.VALUE.value, cls.UNIT.value]


class ComputedSelectorColumns(StrEnum):
    EXPERIMENT = "experiment"
    REMOVE_OUTLIERS = "remove_outliers"
