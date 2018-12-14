<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="html" version="3.2" encoding="iso-8859-1" indent="yes"/>

<xsl:template match="TestResult">
  <html>
    <head>
      <title>AGILE test result</title>
      <link rel="stylesheet" type="text/css" href="unittestreport.css" />
    </head>
    <body>
      <a href="index.html">Back</a><p />
      <xsl:apply-templates select="BuildInfo"/>
      <xsl:apply-templates select="TestSuite"/>
    </body>
  </html>
</xsl:template>

<xsl:template match="TestSuite">
  <xsl:variable name="result" select="normalize-space(@result)" />
  <xsl:choose>
  <xsl:when test="starts-with($result,'p')">
    <h1 class="pass">Testsuite 
    <xsl:value-of select="@result" />:
    <xsl:value-of select="@name" /></h1>
  </xsl:when>
  <xsl:otherwise>
    <h1 class="fail">Testsuite 
    <xsl:value-of select="@result" />:
    <xsl:value-of select="@name" /></h1>
  </xsl:otherwise>
  </xsl:choose>
  <table>
    <tr>
      <th>Assertions passed</th>
      <th>Assertions failed</th>
      <th>Expected failures</th>
      <th>Testcases passed</th>
      <th>Testcases failed</th>
      <th>Testcases skipped</th>
      <th>Testcases aborted</th>
    </tr>
    <tr>
      <td><xsl:value-of select="@assertions_passed" /></td>
      <td><xsl:value-of select="@assertions_failed" /></td>
      <td><xsl:value-of select="@expected_failures" /></td>
      <td><xsl:value-of select="@test_cases_passed" /></td>
      <td><xsl:value-of select="@test_cases_failed" /></td>
      <td><xsl:value-of select="@test_cases_skipped" /></td>
      <td><xsl:value-of select="@test_cases_aborted" /></td>
    </tr>
  </table><p />
  <xsl:apply-templates select="TestSuite"/>
  <table><xsl:apply-templates select="TestCase"/></table>
</xsl:template>

<xsl:template match="TestCase">
  <xsl:if test="position()=1">
    <tr>
      <th class="left">Testcase</th>
      <th>Result</th>
      <th>Assertions passed</th>
      <th>Assertions failed</th>
      <th>Expected failures</th>
    </tr>
  </xsl:if>
  <xsl:variable name="result" select="normalize-space(@result)" />
    <tr>
      <td class="left"><xsl:value-of select="@name" /></td>
  <xsl:choose>
  <xsl:when test="starts-with($result,'p')">
      <td class="pass"><xsl:value-of select="@result" /></td>
  </xsl:when>
  <xsl:otherwise>
      <td class="fail"><xsl:value-of select="@result" /></td>
  </xsl:otherwise>
  </xsl:choose>
      <td><xsl:value-of select="@assertions_passed" /></td>
      <td><xsl:value-of select="@assertions_failed" /></td>
      <td><xsl:value-of select="@expected_failures" /></td>
  </tr>
</xsl:template>

<xsl:template match="BuildInfo">
  <h1> Build Info</h1>
  <table>
    <tr>
      <th>Platform</th>
      <th>Compiler</th>
      <th>STL</th>
      <th>Boost</th>
      <th>Date</th>
      <th>Time</th>
    </tr>
    <tr>
      <td><xsl:value-of select="@platform" /></td>
      <td><xsl:value-of select="@compiler" /></td>
      <td><xsl:value-of select="@stl" /></td>
      <td><xsl:value-of select="@boost" /></td>
      <td><xsl:value-of select="@date" /></td>
      <td><xsl:value-of select="@time" /></td>
    </tr>
  </table>
</xsl:template>
</xsl:stylesheet>
