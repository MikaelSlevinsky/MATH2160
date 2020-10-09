### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 8ec37a68-a8ef-11ea-11cf-ad17550172a4
using Images, LinearAlgebra, Plots

# ╔═╡ 0961cf3a-a8d3-11ea-0d49-2bb1ea933f04
p = load(joinpath(@__DIR__, "Estelada.png"))
#p = load(joinpath(@__DIR__, "MapleLeaf.png"))
#p = load(joinpath(@__DIR__, "StarsAndStripes.gif"))
#p = load(joinpath(@__DIR__, "Communism.jpg"))
#p = process_raw_camera_data(raw_camera_data)

# ╔═╡ 914c04d8-a8d3-11ea-0e25-8bb1c5c224b6
typeof(p)

# ╔═╡ 9741e2ae-a8d3-11ea-2f1e-5bf3985acef9
size(p)

# ╔═╡ 8a202072-a8d3-11ea-1741-dfdf54d7e07d
begin
	m,n = size(p);
	R = Float64[p[i,j].r for i in 1:m, j in 1:n]
	G = Float64[p[i,j].g for i in 1:m, j in 1:n]
	B = Float64[p[i,j].b for i in 1:m, j in 1:n]
end;

# ╔═╡ 87ed6c9c-a8d3-11ea-0172-d342b3bb415f
R

# ╔═╡ 83dfdf4a-a8d3-11ea-0487-311e356935eb
begin
	UR,ΣR,VR = svd(R);
	UG,ΣG,VG = svd(G);
	UB,ΣB,VB = svd(B);
end

# ╔═╡ 7ecd1838-a8d3-11ea-2d13-6f613382e2fc
begin
	σ = ΣR[1] + ΣG[1] + ΣB[1]
	findfirst(i -> i < 1e-3*σ, ΣR),findfirst(i -> i < 1e-3*σ, ΣG),findfirst(i -> i < 1e-3*σ, ΣB)
end

# ╔═╡ 7b72e9e4-a8d3-11ea-0574-49ec485ef91e
begin
	scatter(1:length(ΣR), ΣR; xscale=:log10, yscale=:log10, color=:red, label="Red")
	scatter!(1:length(ΣG), ΣG; color=:green, label="Green")
	scatter!(1:length(ΣB), ΣB; color=:blue, label="Blue")
	xlabel!("Index of singular values")
	ylabel!("Magnitude of singular values")
end

# ╔═╡ f293b524-a8d3-11ea-3c60-017aa28eb3ae
md"""
`r = ` $(@bind r html"<input type='range' max='25'>")
"""

# ╔═╡ 19d55854-a8d4-11ea-014f-01b30e338ba5
begin
	RR = UR[:,1:r]*(ΣR[1:r].*VR[:,1:r]')
    GG = UG[:,1:r]*(ΣG[1:r].*VG[:,1:r]')
    BB = UB[:,1:r]*(ΣB[1:r].*VB[:,1:r]')
    p[:] = [RGB{Normed{UInt8,8}}(clamp(RR[i,j], 0.f0, 1.f0), clamp(GG[i,j], 0.f0, 1.f0), clamp(BB[i,j], 0.f0, 1.f0)) for i in 1:m, j in 1:n]
end

# ╔═╡ 206b2478-a8d4-11ea-051b-4fbabb9838fa
md"""
`rR = ` $(@bind rR html"<input type='range' max='25'>")
`rG = ` $(@bind rG html"<input type='range' max='25'>")
`rB = ` $(@bind rB html"<input type='range' max='25'>")
"""

# ╔═╡ ed11b038-a8d3-11ea-0bae-4f262ef72a82
begin
	RRR = UR[:,1:rR]*(ΣR[1:rR].*VR[:,1:rR]')
    GGG = UG[:,1:rG]*(ΣG[1:rG].*VG[:,1:rG]')
    BBB = UB[:,1:rB]*(ΣB[1:rB].*VB[:,1:rB]')
    p[:] = [RGB{Normed{UInt8,8}}(clamp(RRR[i,j], 0.f0, 1.f0), clamp(GGG[i,j], 0.f0, 1.f0), clamp(BBB[i,j], 0.f0, 1.f0)) for i in 1:m, j in 1:n]
end

# ╔═╡ 2295dd16-f88d-11ea-08b8-4be69d358ee2
function camera_input(;max_size=500, default_url="https://i.imgur.com/SUmi94P.png")
"""
<span class="pl-image waiting-for-permission">
<style>
	
	.pl-image.popped-out {
		position: fixed;
		top: 0;
		right: 0;
		z-index: 5;
	}
	.pl-image #video-container {
		width: 500px;
	}
	.pl-image video {
		border-radius: 1rem 1rem 0 0;
	}
	.pl-image.waiting-for-permission #video-container {
		display: none;
	}
	.pl-image #prompt {
		display: none;
	}
	.pl-image.waiting-for-permission #prompt {
		width: 500px;
		height: 400px;
		display: grid;
		place-items: center;
		font-family: monospace;
		font-weight: bold;
		text-decoration: underline;
		cursor: pointer;
		border: 5px dashed rgba(0,0,0,.5);
	}
	.pl-image video {
		display: block;
	}
	.pl-image .bar {
		width: inherit;
		display: flex;
		z-index: 6;
	}
	.pl-image .bar#top {
		position: absolute;
		flex-direction: column;
	}
	
	.pl-image .bar#bottom {
		background: black;
		border-radius: 0 0 1rem 1rem;
	}
	.pl-image .bar button {
		flex: 0 0 auto;
		background: rgba(255,255,255,.8);
		border: none;
		width: 2rem;
		height: 2rem;
		border-radius: 100%;
		cursor: pointer;
		z-index: 7;
	}
	.pl-image .bar button#shutter {
		width: 3rem;
		height: 3rem;
		margin: -1.5rem auto .2rem auto;
	}
	.pl-image video.takepicture {
		animation: pictureflash 200ms linear;
	}
	@keyframes pictureflash {
		0% {
			filter: grayscale(1.0) contrast(2.0);
		}
		100% {
			filter: grayscale(0.0) contrast(1.0);
		}
	}
</style>
	<div id="video-container">
		<div id="top" class="bar">
			<button id="stop" title="Stop video">✖</button>
			<button id="pop-out" title="Pop out/pop in">⏏</button>
		</div>
		<video playsinline autoplay></video>
		<div id="bottom" class="bar">
		<button id="shutter" title="Click to take a picture">📷</button>
		</div>
	</div>
		
	<div id="prompt">
		<span>
		Enable webcam
		</span>
	</div>
<script>
	// based on https://github.com/fonsp/printi-static (by the same author)
	const span = this.currentScript.parentElement
	const video = span.querySelector("video")
	const popout = span.querySelector("button#pop-out")
	const stop = span.querySelector("button#stop")
	const shutter = span.querySelector("button#shutter")
	const prompt = span.querySelector(".pl-image #prompt")
	const maxsize = $(max_size)
	const send_source = (source, src_width, src_height) => {
		const scale = Math.min(1.0, maxsize / src_width, maxsize / src_height)
		const width = Math.floor(src_width * scale)
		const height = Math.floor(src_height * scale)
		const canvas = html`<canvas width=\${width} height=\${height}>`
		const ctx = canvas.getContext("2d")
		ctx.drawImage(source, 0, 0, width, height)
		span.value = {
			width: width,
			height: height,
			data: ctx.getImageData(0, 0, width, height).data,
		}
		span.dispatchEvent(new CustomEvent("input"))
	}
	
	const clear_camera = () => {
		window.stream.getTracks().forEach(s => s.stop());
		video.srcObject = null;
		span.classList.add("waiting-for-permission");
	}
	prompt.onclick = () => {
		navigator.mediaDevices.getUserMedia({
			audio: false,
			video: {
				facingMode: "environment",
			},
		}).then(function(stream) {
			stream.onend = console.log
			window.stream = stream
			video.srcObject = stream
			window.cameraConnected = true
			video.controls = false
			video.play()
			video.controls = false
			span.classList.remove("waiting-for-permission");
		}).catch(function(error) {
			console.log(error)
		});
	}
	stop.onclick = () => {
		clear_camera()
	}
	popout.onclick = () => {
		span.classList.toggle("popped-out")
	}
	shutter.onclick = () => {
		const cl = video.classList
		cl.remove("takepicture")
		void video.offsetHeight
		cl.add("takepicture")
		video.play()
		video.controls = false
		console.log(video)
		send_source(video, video.videoWidth, video.videoHeight)
	}
	
	
	document.addEventListener("visibilitychange", () => {
		if (document.visibilityState != "visible") {
			clear_camera()
		}
	})
	// Set a default image
	const img = html`<img crossOrigin="anonymous">`
	img.onload = () => {
	console.log("helloo")
		send_source(img, img.width, img.height)
	}
	img.src = "$(default_url)"
	console.log(img)
</script>
</span>
""" |> HTML
end

# ╔═╡ 47f8bf3c-f88e-11ea-36a1-8d9ed51b27eb
@bind raw_camera_data camera_input()

# ╔═╡ 375a2662-f88d-11ea-0127-89f32aeae63f
function process_raw_camera_data(raw_camera_data)
	# the raw image data is a long byte array, we need to transform it into something
	# more "Julian" - something with more _structure_.
	
	# The encoding of the raw byte stream is:
	# every 4 bytes is a single pixel
	# every pixel has 4 values: Red, Green, Blue, Alpha
	# (we ignore alpha for this notebook)
	
	# So to get the red values for each pixel, we take every 4th value, starting at 
	# the 1st:
	reds_flat = UInt8.(raw_camera_data["data"][1:4:end])
	greens_flat = UInt8.(raw_camera_data["data"][2:4:end])
	blues_flat = UInt8.(raw_camera_data["data"][3:4:end])
	
	# but these are still 1-dimensional arrays, nicknamed 'flat' arrays
	# We will 'reshape' this into 2D arrays:
	
	width = raw_camera_data["width"]
	height = raw_camera_data["height"]
	
	# shuffle and flip to get it in the right shape
	reds = reshape(reds_flat, (width, height))' / 255.0
	greens = reshape(greens_flat, (width, height))' / 255.0
	blues = reshape(blues_flat, (width, height))' / 255.0
	
	# we have our 2D array for each color
	# Let's create a single 2D array, where each value contains the R, G and B value of 
	# that pixel
	
	RGB.(reds, greens, blues)
end

# ╔═╡ Cell order:
# ╠═8ec37a68-a8ef-11ea-11cf-ad17550172a4
# ╠═47f8bf3c-f88e-11ea-36a1-8d9ed51b27eb
# ╠═0961cf3a-a8d3-11ea-0d49-2bb1ea933f04
# ╠═914c04d8-a8d3-11ea-0e25-8bb1c5c224b6
# ╠═9741e2ae-a8d3-11ea-2f1e-5bf3985acef9
# ╠═8a202072-a8d3-11ea-1741-dfdf54d7e07d
# ╠═87ed6c9c-a8d3-11ea-0172-d342b3bb415f
# ╠═83dfdf4a-a8d3-11ea-0487-311e356935eb
# ╠═7ecd1838-a8d3-11ea-2d13-6f613382e2fc
# ╠═7b72e9e4-a8d3-11ea-0574-49ec485ef91e
# ╟─f293b524-a8d3-11ea-3c60-017aa28eb3ae
# ╠═19d55854-a8d4-11ea-014f-01b30e338ba5
# ╟─206b2478-a8d4-11ea-051b-4fbabb9838fa
# ╠═ed11b038-a8d3-11ea-0bae-4f262ef72a82
# ╟─2295dd16-f88d-11ea-08b8-4be69d358ee2
# ╟─375a2662-f88d-11ea-0127-89f32aeae63f
