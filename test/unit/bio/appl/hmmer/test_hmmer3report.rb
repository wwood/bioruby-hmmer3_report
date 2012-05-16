#
#
#
# Copyright::   Copyright (C) 2011
#               Christian Zmasek <cmzmasek@yahoo.com>
# License::    The Ruby License
#
#  $Id:$
#
#
# loading helper routine for testing bioruby
require 'helper'

module Bio
  module Testreport

    HMMER_TEST_DATA = Pathname.new(File.join('test','data','HMMER')).cleanpath.to_s

    def self.hmmsearch_domtblout_empty
      File.join HMMER_TEST_DATA, 'hmmsearch_domtblout_empty.out'
    end

    def self.hmmsearch_domtblout
      File.join HMMER_TEST_DATA, 'hmmsearch_domtblout.out'
    end

    def self.hmmsearch_tblout
      File.join HMMER_TEST_DATA, 'hmmsearch_tblout.out'
    end

    def self.hmmscan_domtblout
      File.join HMMER_TEST_DATA, 'hmmscan_domtblout.out'
    end

    def self.hmmscan_tblout
      File.join HMMER_TEST_DATA, 'hmmscan_tblout.out'
    end

  end # Testreport

  class Testreport_class_methods < Test::Unit::TestCase


    def test_hmmsearch_domtblout_empty
      filename = Testreport.hmmsearch_domtblout_empty

      assert_instance_of(Bio::Hmmer3Report,
        report = Bio::Hmmer3Report.new(filename))

      assert_instance_of(Array,
        report.hits)

      assert_equal(0, report.hits.length)
    end #


    def test_hmmsearch_domtblout
      filename = Testreport.hmmsearch_domtblout

      assert_instance_of(Bio::Hmmer3Report,
        report = Bio::Hmmer3Report.new(filename))

      assert_instance_of(Array,
        report.hits)

      assert_instance_of(Bio::HmmerPerDomainHit,
        report.hits[0])

      assert_kind_of(Bio::Hmmer3Hit,
        report.hits[0])

      assert_equal(2,report.hits.length)
      assert_equal("wnt", report.hits[0].target_name)
      assert_equal("PF00110.12", report.hits[0].target_accession)
      assert_equal(310, report.hits[0].target_length)
      assert_equal("sp|P56705|WNT4_HUMAN", report.hits[0].query_name)
      assert_equal("-", report.hits[0].query_accession)
      assert_equal(351, report.hits[0].query_length)
      assert_equal(1.2e-125, report.hits[0].full_sequence_e_value)
      assert_equal(418.5, report.hits[0].full_sequence_score)
      assert_equal(23.4, report.hits[0].full_sequence_bias)
      assert_equal(1, report.hits[0].domain_number)
      assert_equal(1, report.hits[0].domain_sum)
      assert_equal(2.2e-129, report.hits[0].domain_c_e_value)
      assert_equal(1.3e-125, report.hits[0].domain_i_e_value)
      assert_equal(418.3, report.hits[0].domain_score)
      assert_equal(16.2, report.hits[0].domain_bias)
      assert_equal(2, report.hits[0].hmm_coord_from)
      assert_equal(310, report.hits[0].hmm_coord_to)
      assert_equal(43, report.hits[0].ali_coord_from)
      assert_equal(351, report.hits[0].ali_coord_to)
      assert_equal(42, report.hits[0].env_coord_from)
      assert_equal(351, report.hits[0].env_coord_to)
      assert_equal(0.97, report.hits[0].acc)
      assert_equal("wnt family", report.hits[0].target_description)
      assert_equal("MEF2_binding", report.hits[1].target_name)
      assert_equal(0.17, report.hits[1].full_sequence_e_value)
      assert_equal(6.4e-05, report.hits[1].domain_c_e_value)
      assert_equal("MEF2 binding", report.hits[1].target_description)

    end # test_hmmsearch_domtblout

    def test_hmmsearch_tblout

      filename = Testreport.hmmsearch_tblout

      assert_instance_of(Bio::Hmmer3Report,
        report = Bio::Hmmer3Report.new(filename))

      assert_instance_of(Array,
        report.hits)

      assert_instance_of(Bio::HmmerPerSequenceHit,
        report.hits[0])

      assert_kind_of(Bio::Hmmer3Hit,
        report.hits[0])

      assert_equal(1,report.hits.length)
      assert_equal("wnt", report.hits[0].target_name)
      assert_equal("PF00110.12", report.hits[0].target_accession)
      assert_equal("sp|P56705|WNT4_HUMAN", report.hits[0].query_name)
      assert_equal("-", report.hits[0].query_accession)
      assert_equal(9.7e-130, report.hits[0].full_sequence_e_value)
      assert_equal(418.5, report.hits[0].full_sequence_score)
      assert_equal(23.4, report.hits[0].full_sequence_bias)
      assert_equal(1.1e-129, report.hits[0].best_1_domain_e_value)
      assert_equal(418.3, report.hits[0].best_1_domain_score)
      assert_equal(16.2, report.hits[0].best_1_domain_bias)
      assert_equal(1.0, report.hits[0].domain_number_est_exp)
      assert_equal(1, report.hits[0].domain_number_est_reg)
      assert_equal(2, report.hits[0].domain_number_est_clu)
      assert_equal(3, report.hits[0].domain_number_est_ov)
      assert_equal(4, report.hits[0].domain_number_est_env)
      assert_equal(5, report.hits[0].domain_number_est_dom)
      assert_equal(6, report.hits[0].domain_number_est_rep)
      assert_equal(7, report.hits[0].domain_number_est_inc)
      assert_equal("wnt family", report.hits[0].target_description)

    end # test_hmmsearch_tblout

    def test_hmmscan_domtblout

      filename = Testreport.hmmscan_domtblout

      assert_instance_of(Bio::Hmmer3Report,
        report = Bio::Hmmer3Report.new(filename))

      assert_instance_of(Array,
        report.hits)

      assert_instance_of(Bio::HmmerPerDomainHit,
        report.hits[0])

      assert_kind_of(Bio::Hmmer3Hit,
        report.hits[0])

      assert_equal(5,report.hits.length)
      assert_equal("Ion_trans", report.hits[0].target_name)
      assert_equal("PF00520.24", report.hits[0].target_accession)
      assert_equal(201, report.hits[0].target_length)
      assert_equal("jgi|Nemve1|7|gw.28.1.1", report.hits[0].query_name)
      assert_equal("-", report.hits[0].query_accession)
      assert_equal(1604, report.hits[0].query_length)
      assert_equal(6.3e-169, report.hits[0].full_sequence_e_value)
      assert_equal(557.4, report.hits[0].full_sequence_score)
      assert_equal(95.3, report.hits[0].full_sequence_bias)
      assert_equal(1, report.hits[0].domain_number)
      assert_equal(4, report.hits[0].domain_sum)
      assert_equal(1.5e-41, report.hits[0].domain_c_e_value)
      assert_equal(3e-38, report.hits[0].domain_i_e_value)
      assert_equal(130.8, report.hits[0].domain_score)
      assert_equal(11.1, report.hits[0].domain_bias)
      assert_equal(3, report.hits[0].hmm_coord_from)
      assert_equal(171, report.hits[0].hmm_coord_to)
      assert_equal(140, report.hits[0].ali_coord_from)
      assert_equal(307, report.hits[0].ali_coord_to)
      assert_equal(139, report.hits[0].env_coord_from)
      assert_equal(346, report.hits[0].env_coord_to)
      assert_equal(0.81, report.hits[0].acc)
      assert_equal("Ion transport protein", report.hits[0].target_description)
    end # test_hmmscan_domtblout

    def test_hmmscan_tblout
      filename = Testreport.hmmscan_tblout

      assert_instance_of(Bio::Hmmer3Report,
        report = Bio::Hmmer3Report.new(filename))

      assert_instance_of(Array,
        report.hits)

      assert_instance_of(Bio::HmmerPerSequenceHit,
        report.hits[0])

      assert_kind_of(Bio::Hmmer3Hit,
        report.hits[0])

      assert_equal(1,report.hits.length)
      assert_equal("wnt", report.hits[0].target_name)
      assert_equal("PF00110.12", report.hits[0].target_accession)
      assert_equal("sp|P56705|WNT4_HUMAN", report.hits[0].query_name)
      assert_equal("-", report.hits[0].query_accession)
      assert_equal(9.7e-130, report.hits[0].full_sequence_e_value)
      assert_equal(418.5, report.hits[0].full_sequence_score)
      assert_equal(23.4, report.hits[0].full_sequence_bias)
      assert_equal(1.1e-129, report.hits[0].best_1_domain_e_value)
      assert_equal(418.3, report.hits[0].best_1_domain_score)
      assert_equal(16.2, report.hits[0].best_1_domain_bias)
      assert_equal(1.0, report.hits[0].domain_number_est_exp)
      assert_equal(1, report.hits[0].domain_number_est_reg)
      assert_equal(0, report.hits[0].domain_number_est_clu)
      assert_equal(0, report.hits[0].domain_number_est_ov)
      assert_equal(1, report.hits[0].domain_number_est_env)
      assert_equal(1, report.hits[0].domain_number_est_dom)
      assert_equal(1, report.hits[0].domain_number_est_rep)
      assert_equal(1, report.hits[0].domain_number_est_inc)
      assert_equal("wnt family", report.hits[0].target_description)
    end # test_hmmscan_tblout

  end # Testreport _class_methods

end #  Bio
