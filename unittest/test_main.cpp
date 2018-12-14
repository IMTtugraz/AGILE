#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE AGILE tests
#define BOOST_TEST_MAIN 
 
#include <boost/test/unit_test.hpp>

#include <boost/test/output/xml_report_formatter.hpp>
#include <boost/config.hpp>
#include <boost/version.hpp>

class AGILEXmlReportFormatter
  : public boost::unit_test::output::xml_report_formatter
{
  public:
    void results_report_start( std::ostream& ostr)
    {
      // attache the heading for the XML file
      ostr << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>";
      ostr << "<?xml-stylesheet type=\"text/xsl\" "
              "href=\"./unittestreport.xsl\" ?>";
      boost::unit_test::output::xml_report_formatter::results_report_start(
        ostr);

      // some information about this system
      // @see boost/test/impl/xml_log_formatter.ipp
      ostr  << "<BuildInfo"
            << " platform=\"" << BOOST_PLATFORM << '\"'
            << " compiler=\"" << BOOST_COMPILER << '\"'
            << " stl=\""      << BOOST_STDLIB << '\"'
            << " boost=\""    << BOOST_VERSION/100000     << "."
                              << BOOST_VERSION/100 % 1000 << "."
                              << BOOST_VERSION % 100      << '\"'
            << " date=\""     << __DATE__ << '\"'
            << " time=\""     << __TIME__ << '\"'
            << "/>";
    }
};

// Global fixture to set our formatter
struct GlobalConfigurationFixture
{
  GlobalConfigurationFixture()
  {
    boost::unit_test::results_reporter::set_format(
      new AGILEXmlReportFormatter());
  }

  ~GlobalConfigurationFixture() { }
};

BOOST_GLOBAL_FIXTURE( GlobalConfigurationFixture );
