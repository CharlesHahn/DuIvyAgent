"""
DuIvyAgent

TODO
1. 错误选组什么的导致run_terminal不能正常结束，p.wait
2. 增加文档理解能力
3. IMPORTANT 增加任务拆分的逻辑能力
4. 有时候会出现知识库的重复query,不明白为什么？


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
        self.messages = [{"role": "system", "content": system_prompt}]

        self.knowledge = KnowledgeBase(self.config["knowledge_json"])
        all_queries = self.knowledge.get_all_queries()
        print("="*80)
        print("Available queries of knowledge base:")
        for i, query in enumerate(all_queries):
            print(f" [{i+1:>2}] {query}")
        print("="*80)


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
                        },
                        "required": ["command"],
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
            # "max_tokens": self.config["max_tokens"],
            "tools": tools,
        }

        user_prompt = input("\n" + ">"*80 + "\nUser >>> ")
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
            print("*"*80)
            print(result)
            print("*"*80)
            resp_content = result["choices"][0]["message"]["content"]
            if resp_content != "":
                print(f"\nNEO ->>> {resp_content}")
            
            # print(result["choices"][0]["message"])
            
            tool_calls = ""
            tool_calls_cmd = ""
            tool_calls_query = ""
            if "tool_calls" in result["choices"][0]["message"]:
                tool_calls = result["choices"][0]["message"]["tool_calls"][0]["function"]
                tool_calls_args = json.loads(
                    result["choices"][0]["message"]["tool_calls"][0]["function"]["arguments"]
                )
                if tool_calls.get("name", "") == "run_command":
                    tool_calls_cmd = tool_calls_args.get("command", "")
                elif tool_calls.get("name", "") == "query_knowledge":
                    tool_calls_query = tool_calls_args.get("query", "")
                else:
                    pass
            else:
                pass


            self.messages.append({"role": "user", "content": resp_content})
            # self.messages.append(result["choices"][0]["message"])
            # print(self.messages)

            ## parse command
            if tool_calls_cmd != "":
                print(f"NEO >>> Running command: {tool_calls_cmd} (Y/n) ", end="")
                if input().strip() in ["y", ""]:
                    try:
                        cmd_res = run_terminal(tool_calls_cmd)
                        print(f"returncode: {cmd_res['returncode']}")
                        print(f"output: >>> \n {cmd_res['output']}")
                        print(f"error: >>> \n {cmd_res['error']}")
                    except KeyboardInterrupt:
                        print("Command interrupted by user. Try to fix the problem by NEO")
                        cmd_res = "Command interrupted by user. Something wrong with the command, lack of input? please check the command and try again."
                    user_prompt = str(cmd_res)
                else:
                    user_prompt = input("\n"+">"*80 +"\nUser >>> ")
                    if user_prompt == "exit":
                        break
            elif tool_calls_query != "":
                print(f"NEO >>> querying knowledge base for {tool_calls_query}")
                user_prompt = f"The answer to your query is below: \n{'>'*80}\n{str(self.knowledge.query_pair(tool_calls_query))} \n {'>'*80}"
            else:
                user_prompt = input("\n"+">"*80+"\nUser >>> ")
                if user_prompt == "exit":
                    break
            



def main():
    NEO()()


if __name__ == "__main__":
    main()
