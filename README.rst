A Proposal Tutorial for Exoplanet Atmospheres
=============================================

Created for the Sagan Summer School

Setup Up 
--------

Clone the repository 

.. code-block:: bash 

	git clone https://github.com/natashabatalha/propexo.git
	cd propexo

Create conda environment from `yml` file. 

.. code-block:: bash 

	conda env create -f propexo.yml

Table of Contents
-----------------

.. raw:: html

    <embed>
	<div class="toc"><ul class="toc-item"><li><span><a href="#What-accuracy-in-planet-properties-do-I-need?" data-toc-modified-id="What-accuracy-in-planet-properties-do-I-need?-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>What accuracy in planet properties do I need?</a></span><ul class="toc-item"><li><span><a href="#Using-Exoplanet-Archive-API-to-Query-Confirmed-Targets" data-toc-modified-id="Using-Exoplanet-Archive-API-to-Query-Confirmed-Targets-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>Using <code>Exoplanet Archive API</code> to Query Confirmed Targets</a></span></li></ul></li><li><span><a href="#What-tools-do-I-need?-How-do-I-use-them?" data-toc-modified-id="What-tools-do-I-need?-How-do-I-use-them?-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>What tools do I need? How do I use them?</a></span><ul class="toc-item"><li><span><a href="#Use-PandExo-to-run-initial-constant-$(R_p/R_*)^2$-to-determine-approx-precision" data-toc-modified-id="Use-PandExo-to-run-initial-constant-$(R_p/R_*)^2$-to-determine-approx-precision-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>Use <code>PandExo</code> to run initial constant $(R_p/R_*)^2$ to determine approx precision</a></span></li><li><span><a href="#Use-CHIMERA-to-determine-first-guess-atmospheric-transmission-signal" data-toc-modified-id="Use-CHIMERA-to-determine-first-guess-atmospheric-transmission-signal-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>Use <code>CHIMERA</code> to determine first guess atmospheric transmission signal</a></span></li><li><span><a href="#Check-Exo.MAST-for-available-data-so-we-can-validate-our-assumptions" data-toc-modified-id="Check-Exo.MAST-for-available-data-so-we-can-validate-our-assumptions-2.3"><span class="toc-item-num">2.3&nbsp;&nbsp;</span>Check <code>Exo.MAST</code> for available data so we can validate our assumptions</a></span></li><li><span><a href="#Use-PICASO-to-determine-first-guess-atmospheric-emission-signal" data-toc-modified-id="Use-PICASO-to-determine-first-guess-atmospheric-emission-signal-2.4"><span class="toc-item-num">2.4&nbsp;&nbsp;</span>Use <code>PICASO</code> to determine first guess atmospheric emission signal</a></span></li></ul></li><li><span><a href="#How-can-I-&quot;prove&quot;-observability?" data-toc-modified-id="How-can-I-&quot;prove&quot;-observability?-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>How can I "prove" observability?</a></span><ul class="toc-item"><li><span><a href="#Can-an-atmosphere-be-detected:-Addressing-cloud-concerns-and-quantifying-statistical-significance-in-transmission" data-toc-modified-id="Can-an-atmosphere-be-detected:-Addressing-cloud-concerns-and-quantifying-statistical-significance-in-transmission-3.1"><span class="toc-item-num">3.1&nbsp;&nbsp;</span>Can an atmosphere be detected: Addressing cloud concerns and quantifying statistical significance in transmission</a></span></li><li><span><a href="#Can-an-atmosphere-be-detected:-Addressing-unknown-climate-and-quantifying-statistical-significance-in-emission" data-toc-modified-id="Can-an-atmosphere-be-detected:-Addressing-unknown-climate-and-quantifying-statistical-significance-in-emission-3.2"><span class="toc-item-num">3.2&nbsp;&nbsp;</span>Can an atmosphere be detected: Addressing unknown climate and quantifying statistical significance in emission</a></span></li><li><span><a href="#Can-a-specific-molecule-be-detected?" data-toc-modified-id="Can-a-specific-molecule-be-detected?-3.3"><span class="toc-item-num">3.3&nbsp;&nbsp;</span>Can a specific molecule be detected?</a></span></li><li><span><a href="#Can-any-physical-parameters-be-constrained?-Information-content-theory-for-initial-constraint-estimates" data-toc-modified-id="Can-any-physical-parameters-be-constrained?-Information-content-theory-for-initial-constraint-estimates-3.4"><span class="toc-item-num">3.4&nbsp;&nbsp;</span>Can any physical parameters be constrained? Information content theory for initial constraint estimates</a></span></li></ul></li></ul></div>
	</div>
	</div>
	</div>
	</embed>