<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" styleCategories="Symbology">
 <pipe>
  <rasterrenderer type="singlebandpseudocolor" band="1" opacity="0.8" alphaBand="-1" classificationMin="0" classificationMax="1">
   <rasterTransparency/>
   <rastershader>
    <colorrampshader colorRampType="DISCRETE" classificationMode="2" clip="0" minimumValue="0" maximumValue="1" labelPrecision="2">
     <item alpha="255" value="0.2" label="very low" color="#1a9641"/>
     <item alpha="255" value="0.4" label="low" color="#a6d96a"/>
     <item alpha="255" value="0.55" label="moderate" color="#ffffbf"/>
     <item alpha="255" value="0.7" label="high" color="#fdae61"/>
     <item alpha="255" value="inf" label="very high" color="#d7191c"/>
    </colorrampshader>
   </rastershader>
  </rasterrenderer>
  <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
  <rasterresampler maxOversampling="2"/>
 </pipe>
 <blendMode>0</blendMode>
</qgis>
