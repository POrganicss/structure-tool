from openai import OpenAI

client = OpenAI(
  organization='org-HkyNXdToMEQelNCGc2cW4Juq',
) 

stream = client.chat.completions.create(
    model="gpt-4",
    messages=[{"role": "user", "content": "Say this is a test"}],
    stream=True,
)
for chunk in stream:
    if chunk.choices[0].delta.content is not None:
        print(chunk.choices[0].delta.content, end="")