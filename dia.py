"""
DuIvyAgent

"""

import json

import yaml
import requests

from prompt import system_prompt
from knowledge import KnowledgeBase
from terminal import run_terminal


class NEO(object):

    def __init__(self):
        with open("config.yaml", "r") as fo:
            self.config = yaml.safe_load(fo)
        # print(self.config)

        self.knowledge = KnowledgeBase(self.config["knowledge_json"])
        all_queries = self.knowledge.get_all_queries()
        print(all_queries)
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
            {
                "type": "function",
                "function": {
                    "name": "query_knowledge",
                    "description": "Query knowledge base for detailed information by keywords",
                    "parameters": {
                        "type": "object",
                        "properties": {
                            "query": {
                                "type": "string",
                                "description": "the keywords to query",
                                "enum": self.knowledge.get_all_queries(),
                            },
                        },
                        "required": ["query"],
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
            
            print(result["choices"][0]["message"])

            tool_calls = ""
            tool_calls_cmd = ""
            tool_calls_description = ""
            tool_calls_query = ""
            if "tool_calls" in result["choices"][0]["message"]:
                tool_calls = result["choices"][0]["message"]["tool_calls"][0]["function"]
                print(tool_calls)
                tool_calls_args = json.loads(
                    result["choices"][0]["message"]["tool_calls"][0]["function"]["arguments"]
                )
                if tool_calls.get("name", "") == "run_command":
                    tool_calls_cmd = tool_calls_args.get("command", "")
                    tool_calls_description = tool_calls_args.get("description", "")
                elif tool_calls.get("name", "") == "query_knowledge":
                    tool_calls_query = tool_calls_args.get("query", "")
                    print("xxxxxxx NEO querying knowledge base for ", tool_calls_query)
                else:
                    pass
            else:
                pass

            self.messages.append(
                {
                    "role": "assistant",
                    # "content": f"content:{resp_content};\ncommand:{tool_calls_cmd};\ncommand_description:{tool_calls_description}",
                    "content": str(result["choices"][0]["message"])
                }
            )
            # print(self.messages)

            ## parse command
            if tool_calls_cmd != "":
                print(f"\nNEO >>> {tool_calls_description}")
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
            elif tool_calls_query != "":
                print(f"NEO querying knowledge base for {tool_calls_query}")
                user_prompt = f"The answer to your query is below: \n{'>'*80}\n{str(self.knowledge.query_pair(tool_calls_query))} \n {'>'*80}"
            else:
                user_prompt = input("\n"+">"*78+"\nUser >>> ")
                if user_prompt == "exit":
                    break
            



def main():
    NEO()()


if __name__ == "__main__":
    main()
