<?xml version='1.0' encoding='utf-8'?>

<!-- Stylesheet to upgrade from Mitsuba version 0.4.x to 0.6.0 scenes -->

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">
    <xsl:output method="xml" indent="yes" encoding="utf-8"/>
    <xsl:preserve-space elements="*"/>

    <!-- Update the scene version -->
    <xsl:template match="scene/@version">
        <xsl:attribute name="version">0.6.0</xsl:attribute>
    </xsl:template>

    <!-- Default copy rule -->
    <xsl:template match="@*|node()">
        <xsl:copy>
            <xsl:apply-templates select="@*|node()"/>
        </xsl:copy>
    </xsl:template>
</xsl:stylesheet>
