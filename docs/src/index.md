# GreenFunctionMonteCarlo

```@contents
```

```@eval
using Markdown
function remove_comments(input::String)::String
    s1 = replace(input, r"<!--.*?-->" => "")
end
function strip_to_start(input::String)::String
    start_index = findfirst(r"<!-- START README \(DO NOT DELETE THIS LINE!\) -->", input)
    return start_index === nothing ? "" : input[last(start_index)+1:end]
end

Markdown.parse(remove_comments(strip_to_start(read("../../README.md", String))))
```