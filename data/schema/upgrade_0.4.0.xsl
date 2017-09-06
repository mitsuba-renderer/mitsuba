<?xml version='1.0' encoding='utf-8'?>

<!-- Stylesheet to upgrade from Mitsuba version 0.3.x to 0.4.0 scenes -->

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
    <xsl:output method="xml" indent="yes" encoding="utf-8"/>
    <xsl:preserve-space elements="*"/>

    <!-- Update the scene version -->
    <xsl:template match="scene/@version">
        <xsl:attribute name="version">0.4.0</xsl:attribute>
    </xsl:template>

    <!-- Cameras have been renamed to sensors -->
    <xsl:template match="camera">
        <sensor>
            <xsl:apply-templates select="@*"/>

            <!-- Handle the mapSmallerSide parameter -->
            <xsl:if test="@type='perspective'">
                <string name="fovAxis">
                    <xsl:attribute name="value">
                        <xsl:choose>
                            <xsl:when test="boolean[@name='mapSmallerSide' and @value='false']">
                                <xsl:text>larger</xsl:text>
                            </xsl:when>
                            <xsl:otherwise>
                                <xsl:text>smaller</xsl:text>
                            </xsl:otherwise>
                        </xsl:choose>
                    </xsl:attribute>
                </string>
            </xsl:if>

            <xsl:apply-templates select="node()[not(@name) or @name!='mapSmallerSide']"/>
        </sensor>
    </xsl:template>

    <!-- Fix the handedness yet once more.. -->
    <xsl:template match="camera/transform[@name='toWorld']">
        <transform>
            <xsl:apply-templates select="@*"/>
            <scale x="-1"/>
            <xsl:apply-templates select="node()"/>
        </transform>
    </xsl:template>

    <!-- Rename the 'intensity' parameter of certain luminaires to 'radiance' -->
    <xsl:template match="luminaire[@type='area' or @type='constant']/node()[@name='intensity']/@name">
        <xsl:attribute name="name">radiance</xsl:attribute>
    </xsl:template>

    <!-- Update the 'intensity' parameter of directional light sources -->
    <xsl:template match="luminaire[@type='directional']/node()[@name='intensity']/@name">
        <xsl:attribute name="name">irradiance</xsl:attribute>
    </xsl:template>

    <!-- Rename the luminaireSamples parameter of the direct sampling strategy -->
    <xsl:template match="integrator/node()[@name='luminaireSamples']/@name">
        <xsl:attribute name="name">emitterSamples</xsl:attribute>
    </xsl:template>

    <!-- Rename the depth parameter of samplers -->
    <xsl:template match="sampler/node()[@name='depth']/@name">
        <xsl:attribute name="name">dimension</xsl:attribute>
    </xsl:template>

    <!-- Update the name of the errctrl plugin -->
    <xsl:template match="integrator[@type='errctrl']/@type">
        <xsl:attribute name="type">adaptive</xsl:attribute>
    </xsl:template>

    <!-- Update the name of the exrfilm plugin -->
    <xsl:template match="film[@type='exrfilm']/@type">
        <xsl:attribute name="type">hdrfilm</xsl:attribute>
    </xsl:template>

    <!-- Translate the 'alpha' parameter in the old films -->
    <xsl:template match="film/boolean[@name='alpha']">
        <xsl:choose>
            <xsl:when test="@value='true'">
                <string name="pixelFormat" value="rgba"/>
            </xsl:when>
            <xsl:otherwise>
                <string name="pixelFormat" value="rgb"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>

    <!-- Update the name of the pngfilm plugin -->
    <xsl:template match="film[@type='pngfilm']/@type">
        <xsl:attribute name="type">ldrfilm</xsl:attribute>
    </xsl:template>

    <!-- Update the 'focusDepth' attribute name -->
    <xsl:template match="float[@name='focusDepth']/@name">
        <xsl:attribute name="name">focusDistance</xsl:attribute>
    </xsl:template>

    <!-- Update the 'intensityScale' attribute name -->
    <xsl:template match="float[@name='intensityScale']/@name">
        <xsl:attribute name="name">scale</xsl:attribute>
    </xsl:template>

    <!-- Update the 'densityMultiplier' attribute name -->
    <xsl:template match="float[@name='densityMultiplier']/@name">
        <xsl:attribute name="name">scale</xsl:attribute>
    </xsl:template>

    <!-- Update the blackbody 'multiplier' attribute name -->
    <xsl:template match="blackbody/@multiplier">
        <xsl:attribute name="scale"><xsl:value-of select="."/></xsl:attribute>
    </xsl:template>

    <!-- Luminaires have been renamed to emitters -->
    <xsl:template match="luminaire">
        <emitter>
            <xsl:apply-templates select="@*|node()"/>
        </emitter>
    </xsl:template>

    <!-- Default copy rule -->
    <xsl:template match="@*|node()">
        <xsl:copy>
            <xsl:apply-templates select="@*|node()"/>
        </xsl:copy>
    </xsl:template>
</xsl:stylesheet>
