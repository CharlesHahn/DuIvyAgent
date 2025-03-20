"""
DuIvyAgent

"""

import json

import yaml
import requests

from prompt import system_prompt
from terminal import run_terminal


class NEO(object):

    def __init__(self):
        with open("config.yaml", "r") as fo:
            self.config = yaml.safe_load(fo)
        # print(self.config)
        self.messages = [{"role": "system", "content": system_prompt}]

    def __call__(self):

        tools = [
            {
                "type": "function",
                "function": {
                    "name": "run_command",
                    "description": "Run a shell command",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "command": {
                                "type": "string",
                                "description": "Shell command to run",
                            },
                            "description": {
                                "type": "string",
                                "description": "description of the command",
                            },
                        },
                        "required": ["command", "description"],
                    },
                },
            },
        ]

        headers = {
            "Authorization": f"Bearer {self.config['api_key']}",
            "Content-Type": "application/json",
        }
        data = {
            "model": self.config["model"],
            "temperature": self.config["temperature"],
            "max_tokens": self.config["max_tokens"],
            "tools": tools,
        }

        user_prompt = input("\n>>>>>>>>>\nUser >>> ")
        while True:
            if user_prompt == "exit":
                break

            self.messages.append({"role": "user", "content": user_prompt})
            data["messages"] = self.messages
            # print(self.messages)
            response = requests.post(
                url=self.config["base_url"],
                headers=headers,
                json=data,
                # timeout=self.config["timeout"],
            )
            response.raise_for_status()
            result = response.json()
            resp_content = result["choices"][0]["message"]["content"]
            if resp_content != "":
                print(f"\nNEO >>> {resp_content}")

            if "tool_calls" in result["choices"][0]["message"]:
                tool_calls = json.loads(
                    result["choices"][0]["message"]["tool_calls"][0]["function"][
                        "arguments"
                    ]
                )
                # print(tool_calls)
                tool_calls_cmd = tool_calls.get("command", "")
                tool_calls_description = tool_calls.get("description", "")
                print(f"\nNEO >>> {tool_calls_description}")
            else:
                tool_calls = ""
                tool_calls_cmd = ""
                tool_calls_description = ""
            self.messages.append(
                {
                    "role": "assistant",
                    "content": f"content:{resp_content};\ncommand:{tool_calls_cmd};\ncommand_description:{tool_calls_description}",
                }
            )
            # print(self.messages)

            ## parse command
            if tool_calls_cmd != "":
                print(f"Running command: {tool_calls_cmd} (Y/n) ", end="")
                if input().strip() in ["y", ""]:
                    cmd_res = run_terminal(tool_calls_cmd)
                    print(f"returncode: {cmd_res['returncode']}")
                    print(f"output: >>> \n {cmd_res['output']}")
                    print(f"error: >>> \n {cmd_res['error']}")
                    user_prompt = str(cmd_res)
                else:
                    user_prompt = input("\n"+">"*78 +"\nUser >>> ")
                    if user_prompt == "exit":
                        break
            else:
                user_prompt = input("\n"+">"*78+"\nUser >>> ")
                if user_prompt == "exit":
                    break


def main():
    NEO()()


if __name__ == "__main__":
    main()
