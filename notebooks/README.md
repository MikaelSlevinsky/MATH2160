# Pluto notebooks

Please visit [JuliaLang/Downloads](https://julialang.org/downloads/) to download and install Julia.

**Tip:** platform-specific [instructions](https://julialang.org/downloads/platform/) can help with post-installation conveniences such as launching Julia by typing `julia` in the terminal.

Julia consists of a `Base` library and is also shipped with a *standard set of libraries*. One example of these is `LinearAlgebra`.

To use a standard library package such as `LinearAlgebra`, just type `using LinearAlgebra` at the `julia>` prompt. To see what a package contains, type `LinearAlgebra.<TAB>`. This will show you all the types and functions that are declared within.

To add a [registered package](https://github.com/JuliaRegistries/General) from the Julia community's ecosystem hosted on `GitHub`, switch from execution mode (`julia>`) to package mode (`(@v1.5) pkg> `) by typing `]`. Add [Pluto.jl](https://github.com/fonsp/Pluto.jl) by typing `add Pluto` followed by `enter`. This may take a while as it initializes your package environment for the first time. To exit package mode, press `backspace`. Let's also add `Plots` as it's used in the first notebook.

To use `Pluto`, type `using Pluto; Pluto.run()`. This will open the `Pluto` landing page in your default browser (and by default at the local web address `http://localhost:1234`). Here, you can take a look at a few sample notebooks, start a new notebook, and even open a notebook that's hosted on the Internet by pasting the link in the `Open from file` box.

For a short YouTube video that shows these steps, please see [here](https://m.youtube.com/watch?v=OOjKEgbt8AI).

Our Pluto notebooks can also be viewed *statically* by downloading the PDF file or by using `GitHub`'s [HTML preview](https://htmlpreview.github.io):

- [Introduction.jl](https://htmlpreview.github.io/?https://github.com/MikaelSlevinsky/MATH2160/blob/master/notebooks/Introduction.jl.html)
- [Norms.jl](https://htmlpreview.github.io/?https://github.com/MikaelSlevinsky/MATH2160/blob/master/notebooks/Norms.jl.html)
- [Numbers.jl](https://htmlpreview.github.io/?https://github.com/MikaelSlevinsky/MATH2160/blob/master/notebooks/Numbers.jl.html)
- [Conditioning.jl](https://htmlpreview.github.io/?https://github.com/MikaelSlevinsky/MATH2160/blob/master/notebooks/Conditioning.jl.html)
- [Complexity.jl](https://htmlpreview.github.io/?https://github.com/MikaelSlevinsky/MATH2160/blob/master/notebooks/Complexity.jl.html)
- [RandomPolygons.jl](https://htmlpreview.github.io/?https://github.com/MikaelSlevinsky/MATH2160/blob/master/notebooks/RandomPolygons.jl.html)
- [Svd.jl](https://htmlpreview.github.io/?https://github.com/MikaelSlevinsky/MATH2160/blob/master/notebooks/Svd.jl.html)
- [OrthogonalPolynomials.jl](https://htmlpreview.github.io/?https://github.com/MikaelSlevinsky/MATH2160/blob/master/notebooks/OrthogonalPolynomials.jl.html)
