import sys


class Questions:
    @staticmethod
    def inputEscape(question):
        answer = input(question)
        if answer == "":
            if Questions.askYesorNo("EXIT PROGRAM?"):
                print("EXITING PROGRAM")
                sys.exit()
        return answer

    @staticmethod
    def askMultipleChoice(question, choices):
        if len(choices) == 1:
            return list(choices)[0]
        choice_mapping = {f"{i}": choice for i, choice in enumerate(choices)}
        options = "\n".join([f"{i}: {choice}" for i, choice in choice_mapping.items()])
        choice = Questions.inputEscape(f"{question}\n{options}\n")
        while choice not in choice_mapping.keys():
            choice = Questions.inputEscape(
                f"""Invalid choice, possibilities are:\{options}\n"""
            )
        return choice_mapping[choice]

    @staticmethod
    def askSelectParameter(data, column):
        options = set(data[column])
        print(options)
        answer = Questions.inputEscape(f"""Select {column}?\n{', '.join(options)}\n""")
        while answer not in options:
            answer = Questions.inputEscape(
                f"Invalid choice, possibilities are:\n{', '.join(options)}\n"
            )
        return answer

    @staticmethod
    def askYesorNo(question):
        answer = Questions.inputEscape(f"""{question} (y/n)\n""").upper()
        while answer not in ["Y", "N"]:
            answer = Questions.inputEscape(
                "Invalid choice, possibilities are: (y/n)\n"
            ).upper()
        return answer == "Y"
