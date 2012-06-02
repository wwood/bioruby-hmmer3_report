require 'helper'

module Bio
  class TestDefaultReport < Test::Unit::TestCase
    def test_no_format_specified
      reports = Bio::HMMER::HMMER3.reports(File.open(File.join(HMMER_TEST_DATA, 'test637000001.hmmsearch.txt')))
      assert_kind_of Array, reports
      assert_equal 4, reports.length
      assert_kind_of Bio::HMMER::HMMER3::DefaultHMMSearchReport, reports[0]
    end
    
    def test_tblout_format_specified
      report = Bio::HMMER::HMMER3.reports(File.open(File.join(HMMER_TEST_DATA, 'hmmsearch_tblout.out')), :format => :tblout)
      assert_kind_of Bio::HMMER::HMMER3::TabularReport, report
      assert_equal :tblout, report.format
    end
    
    def test_tblout_format_specified
      report = Bio::HMMER::HMMER3.reports(File.open(File.join(HMMER_TEST_DATA, 'hmmsearch_domtblout.out')), :format => :domtblout)
      assert_kind_of Bio::HMMER::HMMER3::TabularReport, report
      assert_equal :domtblout, report.format
    end
  end
end