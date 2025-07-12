from enum import StrEnum


class DatasetColumn(StrEnum):
    VALUE = "value"
    UNIT = "unit"
    SUBJECT_ID = "subject_id"
    GROUP_ID = "group_id"

    @classmethod
    def get_mandatory_columns(cls):
        return [cls.SUBJECT_ID, cls.VALUE, cls.UNIT]
