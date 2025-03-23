from collections import namedtuple


def create_measurement_class(name, attributes, string_template):

    cls = namedtuple(name, attributes)

    def __eq__(self, other):
        if isinstance(other, cls):
            return other == cls
        if isinstance(other, dict):
            for attribute in other:
                if attribute not in attributes:
                    return False
                if other[attribute] != getattr(self, attribute):
                    return False
            return True
        return False

    cls.__eq__ = __eq__

    def __str__(self):
        string = string_template
        for attribute in attributes:
            string = string.replace(attribute, getattr(self, attribute))
        return string

    cls.__str__ = __str__

    return cls
