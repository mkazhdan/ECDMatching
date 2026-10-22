<center><h2>Extended Convolution Descriptor Matching (Version 1.00)</h2></center>
<center>
<a href="#LINKS">links</a>
<a href="#EXECUTABLES">executables</a>
<a href="#EXAMPLES">examples</a>
<a href="#COMPILATION">compilation</a>
<a href="#CHANGES">changes</a>
</center>
<hr>
This software supports the comparison of two images through the computation and comparison of the dominant Extended Convolution Descriptors.
<hr>
<a name="LINKS"><b>LINKS</b></a><br>
<ul>
<b>Papers:</b>
<a href="https://www.cs.jhu.edu/~misha/MyPapers/ArXiv20.pdf">[Mitchel, Brown, Koller, Weyrich, Rusinkiewicz, and Kazhdan, 2020]</a>
<br>
<b>Executables: </b>
<a href="https://www.cs.jhu.edu/~misha/Code/ECD/Version1.00/ECD.x64.zip">Win64</a><br>
<b>Source Code:</b>
<a href="https://www.cs.jhu.edu/~misha/Code/ECD/Version1.00/ECD.Source.zip">ZIP</a> <a href="https://github.com/mkazhdan/ECD">GitHub</a><br>
<B>Data:</B>
<A HREF="https://www.cs.jhu.edu/~misha/Code/ECD/ECD.Data.zip">ZIP</A><br>
<!--
<b>Older Versions:</b>
<a href="https://www.cs.jhu.edu/~misha/Code/ECD/Version1.10/">V1.00</a>,
<a href="https://www.cs.jhu.edu/~misha/Code/ECD/Version1.00/">V1.00</a>
-->
</ul>
<hr>
<a name="EXECUTABLES"><b>EXECUTABLES</b></a><br>
<ul>

<dl>
<details>
<summary>
<font size="+1"><b>ECDMatching</b></font>:
The executable computes the dominant descriptors for the source and target and uses those to measure the image similarity. Concretely, for each target descriptor it computes the distance to the closest source descriptor and outputs the sum of the distances, over all target descriptors.
</summary>
<dt><b>--in</b> &lt;<i>input source image</i>&gt; &lt;<i>input target image</i>&gt;</dt>
<dd> These two strings specify the the names of the files containing the source and target images.<BR>
The executable supports JPEG and PNG image formats.<BR>
</dd>

<dt>[<b>--mask</b> &lt;<i>input source mask</i>&gt; &lt;<i>input target mask</i>&gt;]</dt>
<dd> These optional string specify masking image files indicating the pixels at which descriptors are to be computed. (Descriptors will not be computed at pixels where the mask file is black.)<BR>
The source (resp. target) mask must have the same resolution as the source (resp. target) image.
</dd>

<dt>[<b>--ecdHalfRes</b> &lt;<i>descriptor half resolution</i>&gt;]</dt>
<dd> This optional integer specifies half the resolution of the ECD descriptor.<BR>
The default value is 5.
</dd>

<dt>[<b>--ecdRadius</b> &lt;<i>descriptor radius (in pixels)</i>&gt;]</dt>
<dd> This optional floating point value specifies radius over which the Extended Convolution Descriptor is computed.<BR>
The default value is twice the ECD half-resolution.
</dd>

<dt>[<b>--ecdNum</b> &lt;<i>source ECD num</i>&gt; &lt;<i>target ECD num</i>&gt;]</dt>
<dd> These two optional integer values specifies the number of dominant descriptors to be computed for the source and target respectively.<BR>
The default values are -1 for the source and 1000 for the target, indicating that the 1000 most dominant target descriptors, should be matched against all the source descriptors to measure target-to-source similarity.
</dd>

<dt>[<b>--stratified</b> &lt;<i>source stratification resolution</i>&gt; &lt;<i>target stratification resolution num</i>&gt;]</dt>
<dd> These two optional integer values specify the stratification resolutions (in pixels) for the source and target. The initial set of descriptors are compute by partitioning the image into squares of the prescribed resolution, and selecting the dominant descriptor within each square.<BR>
A resolution of 0 (or 1) indicates that no stratification should be used.<BR>
The default values are 2 for the source and 8 for the target.
</dd>

</details>
</dl>

</ul>


<hr>
<a name="EXAMPLES"><b>EXAMPLES (WITH SAMPLE DATA)</b></a><br>
For testing purposes, <A HREF="https://www.cs.jhu.edu/~misha/Code/ECD/ECD.Data.zip">the provided archive</A> contains views of the texture-mapped "vase" model from nine different perspectives, including renderings obtained using the reference texture and renderings obtained using synthesized textures.<BR>
The renderings are obtained using the model and photographs generously shared by Paolo Cignoni and Marco Callieri, <A HREF="https://vcg.isti.cnr.it/">Visual Computing Laboratory, ISTI-CNR</A>.

<ul>

<dl>
<details>
<summary>
<font size="+1"><b>ECDMatching</b></font>
</summary>
To run this executable you must specify the source and the target images (as well as associated masks). Assuming you are running from within one of the "view" folders:
<blockquote><code>% ECDMatching --in reference.png 1.png --mask mask.png mask.png</code></blockquote>
This computes the Extended Convolution descriptor for the source and target and gives the sum of the distances from the target descriptors to their closest source descriptors, and prints the sum to the standard output.<BR>
To compare multiple targets to the source, the second argument of the <B>--in</B> (and <B>--mask</B>) parameter can be a ':'-seperated list of target image (and associated mask) file-names:
<blockquote><code>% ECDMatching --in reference.png 1.png:2.png:3.png:4.png --mask mask.png mask.png:mask.png:mask.png:mask.png</code></blockquote>
In this case the sum of distances from each target will be written out to separate lines of the standard output.
</details>
</dl>

</ul>


<hr>
<details>
<summary>
<a name="COMPILATION"><b>COMPILATION</b></a><br>
</summary>
<UL>
<LI>Compilation requires the ZLIB, JPEG, and PNG libraries. These are provided with the Visual Studios solution file. The Makefile assumes that the libraries and associated headers are located in the associated paths used by the compiler.
</UL>
</details>


<hr>
<details>
<summary>
<a name="CHANGES"><b>HISTORY OF CHANGES</b></a><br>
</summary>

<a href="https://www.cs.jhu.edu/~misha/Code/ECDMatching/Version1.00/">Version 1.00</a>:
<ul>
<li> initial source code
</ul>

</details>