require 'helper'

module Bio
  class TestDefaultReport < Test::Unit::TestCase
    def test_splitting
      reports = Bio::HMMER::HMMER3.reports(File.open(File.join(HMMER_TEST_DATA, 'test637000001.hmmsearch.txt')))
      assert_equal 4, reports.length
      assert_kind_of Bio::HMMER::HMMER3::DefaultHMMSearchReport, reports[0]
    end
    
    def test_alignment_when_no_hits
      reports = Bio::HMMER::HMMER3.reports(File.open(File.join(HMMER_TEST_DATA, 'test637000001.hmmsearch.txt')))
      assert_equal [], reports[0].hits
    end

    def test_alignment_when_three_hits
      reports = Bio::HMMER::HMMER3.reports(File.open(File.join(HMMER_TEST_DATA, 'test637000001.hmmsearch.txt')))
      hits = reports[3].hits
      assert_kind_of Array, hits
      
      assert_equal 3, hits.length
      
      
      h = hits[1]
      assert_kind_of Bio::HMMER::HMMER3::DefaultHMMSearchReport::Hit, h
      assert_equal 1, h.hsps.length
      assert_equal '637984252  Acid345_2236 D-isomer specific 2-hydroxyacid dehydrogenase, NAD-binding [Korebacter versatilis Ellin345]',
        h.sequence_name
      
      d = h.hsps[0]
      assert_kind_of Bio::HMMER::HMMER3::DefaultHMMSearchReport::Hit::Hsp, d
      assert_equal 1, d.number
      assert_equal 66.8, d.score
      assert_equal 0.0, d.bias
      assert_equal 2.4e-22, d.c_evalue
      assert_equal 3.8e-19, d.i_evalue
      assert_equal 7, d.hmmfrom
      assert_equal 102, d.hmm_to
      assert_equal 14, d.alifrom
      assert_equal 146, d.ali_to
      assert_equal 8, d.envfrom
      assert_equal 353, d.env_to
      assert_equal 0.77, d.acc
      
      assert_equal 'lreeelellke.gvevevkd...ellteellekakd.adalivrsntkvtaevlekl.pkLkviatagvGvDniDldaakerGIlVtnvpgystesvAE'+
        'la...............................f', d.hmmseq
      assert_equal 'IGKPALERLRAaGYDVEVYPqadPPPKSLIIEKVASgIDGLITTLRDKIDAEVFEAGkGNLKVVAQIAVGFDNINRADANKYKVPFTHTADVLTEATAE'+
        'FAffimaaaarklwtaernvrdlkwgtwhpflpF', d.flatseq
        
      assert_equal 'IVVAEKIAKAAIDLFKQdpTWNVVTPDqVAQKEQLLEQLKGADALIVRSAVFVDAAMLEHADQLRVIGRAGVGVDNIELEAATRKGIAVMNTPGANAIA'+
        'VAEHTiglmlalarfipratetmhagkwekkslqgtelrgktlgivglgriglevarraasfgmtlvahdpyvspaiahdakirladrdevlavadyit'+
        'lhvgltpqtanminattlatmkkgvrivncargeliddaalaeavksghvggaaldvfteeplkaspyhgvpnviltphigGSTAEAQDAVGVQIAHQV'+
        'RDYLQRGVVQNAVN', hits[0].hsps[0].flatseq
        
      assert_equal 'eellekakdadalivrsntkvtaevleklpkLkviatagvGvDniDldaakerGIlVtnvpgystesvAEla...........................'+
        '...................................................................................................'+
        '........................................................fateeaqenmaeeaaenlvaflkgespanav', hits[2].hsps[0].hmmseq
    end
    
    def test_multi_domain_hit_simple
      hits = Bio::HMMER::HMMER3.reports(File.open(File.join(HMMER_TEST_DATA, 'hmmer3multidomainHitSimple.txt')))[0].hits
      assert_equal 1, hits.length
      assert_equal 2, hits[0].hsps.length
      
      assert_equal 'FIGNRIGTFSVLNVIRVMQEMDLSIEDVDALTGSAVGWPkSATFRTIDLVGLDILGHVVGNMKQNVTDErsDLQIPDFYKQMLERKWLGDKTKGGFYK', hits[0].hsps[0].flatseq
      assert_equal 'DTIVEIDAAMRMGFNWEMGPFELWDAAGVEATVGRMKA', hits[0].hsps[1].flatseq
    end
    
    def test_multi_domain_hit
      hits = Bio::HMMER::HMMER3.reports(File.open(File.join(HMMER_TEST_DATA, 'hmmer3multidomainHit.txt')))[0].hits
      assert_equal 4, hits.length
      assert_equal 1, hits[0].hsps.length
      assert_equal 2, hits[2].hsps.length
      
      assert_equal 'FIGNRIGTFSVLNVIRVMQEMDLSIEDVDALTGSAVGWPkSATFRTIDLVGLDILGHVVGNMKQNVTDErsDLQIPDFYKQMLERKWLGDKTKGGFYK', hits[2].hsps[0].flatseq
      assert_equal 'DTIVEIDAAMRMGFNWEMGPFELWDAAGVEATVGRMKA', hits[2].hsps[1].flatseq
    end
  end
end