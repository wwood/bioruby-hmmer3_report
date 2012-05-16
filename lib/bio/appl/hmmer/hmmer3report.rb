#
# = bio/appl/hmmer/hmmer3report.rb - hmmscan/hmmsearch parser
#
# Copyright::   Copyright (C) 2011
#               Christian Zmasek <cmzmasek@yahoo.com>
# License::     The Ruby License
#
# $Id:$
#
# == Description
#
# Parser class for hmmsearch and hmmscan in the HMMER 3 package.
#
# == Examples
#
#    # Input from file:
#    report = Bio::Hmmer3Report.new('/path/to/hmmer_domtblout.output')
#    report.hits.each do |hit|
#      puts hit.target_name
#      puts hit.target_accession
#      puts hit.query_name
#      puts hit.query_accession
#      puts hit.query_length
#      puts hit.full_sequence_e_value
#      puts hit.full_sequence_score
#      puts hit.domain_number
#      puts hit.domain_sum
#      puts hit.domain_c_e_value
#      puts hit.domain_i_e_value
#      puts hit.domain_score
#      puts hit.domain_bias
#      puts hit.hmm_coord_from
#      puts hit.hmm_coord_to
#      puts hit.ali_coord_from
#      puts hit.ali_coord_to
#      puts hit.env_coord_from
#      puts hit.env_coord_to
#      puts hit.acc
#      puts hit.target_description
#    end
#
#
#
#    # Input from string:
#    data = String.new
#    data         << '#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord'
#    data << "\n" << '# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target'
#    data << "\n" << '#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------'
#    data << "\n" << 'Bcl-2                PF00452.13   101 sp|P10415|BCL2_HUMAN -            239   3.7e-30  103.7   0.1   1   1   7.9e-34   4.9e-30  103.3   0.0     1   101    97   195    97   195 0.99 Apoptosis regulator proteins, Bcl-2 family'
#    data << "\n" << 'BH4                  PF02180.11    27 sp|P10415|BCL2_HUMAN -            239   3.9e-15   54.6   0.1   1   1   1.3e-18   8.2e-15   53.6   0.1     2    26     8    32     7    33 0.94 Bcl-2 homology region 4'
#    data << "\n"
#
#    report = Bio::Hmmer3Report.new(data)
#    report.hits.each do |hit|
#      puts hit.target_name
#      puts hit.full_sequence_e_value
#    end
#
#
#
#    # RDF/XML output:
#    puts report.to_rdf
#
#
#
# == References
#
# * HMMER
#   http://hmmer.janelia.org/
#

require 'stringio'

module Bio

  class Hmmer3Report
    def initialize(hmmer_output, format = nil)
      @hits = Array.new
      @line_number = 0
      @format = format
      if File.exists?(hmmer_output.to_s)
        my_hmmer_output_file = File.new(hmmer_output.to_s)
        my_hmmer_output_file.each_line() { |line| parse_line(line) }
      else
        str = StringIO.new(hmmer_output)
        str.each_line() { |line| parse_line(line) }
      end
    end

    attr_reader :hits
    attr_reader :format

    def to_rdf(type = :xml)
      output = String.new
      output << rdf_header(type)
      hits.each do |hit|
        output << hit.to_rdf(type)
      end
      output << rdf_end(type)
      output
    end

    private

    def parse_line(line)
      @line_number += 1
      if line  =~ /^#.+this\s+domain/
        @format = :domtblout
      elsif line =~ /^#.+best\s+1\s+domain/
        @format = :tblout
      elsif line =~ /\S/ && line !~ /^#/
        if @format == nil
          if looks_like_per_domain_result?(line)
            @format = :domtblout
          else
            @format = :tblout
          end
        end
        if @format == :domtblout
          @hits << HmmerPerDomainHit.new(line, @line_number)
        elsif @format == :tblout
          @hits << HmmerPerSequenceHit.new(line, @line_number)
        else
          raise ArgumentError, "attempt to parse hmmscan/hmmsearch output style other than \"domtblout\" or \"tblout\""
        end
      end
    end

    def looks_like_per_domain_result?(line)
      line =~ /^(\S*)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s*(.*)/
    end

    #RDF output
    def rdf_header(type = :xml)
      '<?xml version="1.0"?>' +
        "\n" + '<rdf:RDF' +
        "\n" + ' xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"' +
        "\n" + ' xmlns:' + Hmmer3Hit::RDF_HMMER_HIT+'="http://www.open-bio.org/core/'+Hmmer3Hit::RDF_HMMER_HIT+'/"' +
        "\n" + '>'
    end

    def rdf_end(type = :xml)
      "\n</rdf:RDF>"
    end

  end # class Hmmer3Report

  class Hmmer3Hit
    def initialize
      # This is an abstract class. Prevents 'new' being called on this class
      # and force implementation of 'initialize' in inheriting classes.
      raise NotImplementedError
    end
    attr_reader :target_name
    attr_reader :target_accession
    attr_reader :target_description
    attr_reader :query_name
    attr_reader :query_accession
    attr_reader :full_sequence_e_value
    attr_reader :full_sequence_score
    attr_reader :full_sequence_bias

    # RDF output
    RDF_HMMER_PFAM_RESOURCE = 'http://pfam.sanger.ac.uk/family?acc='
    RDF_HMMER_HIT = 'hmmer_hit'
    RDF_HMMER_TARGET_NAME = 'target_name'
    RDF_HMMER_QUERY_NAME = 'query_name'
    RDF_HMMER_QUERY_ACCESSION = 'query_accession'
    RDF_HMMER_FULL_SEQUENCE_E_VALUE = 'full_sequence_e_value'
    RDF_HMMER_FULL_SEQUENCE_SCORE = 'full_sequence_score'
    RDF_HMMER_FULL_SEQUENCE_BIAS = 'full_sequence_bias'
    RDF_HMMER_TARGET_DESCRIPTION = 'target_description'

    private

    def to_rdfxml( label_a, label_b, value, indent, attribute = '' )
      if attribute && attribute.length > 0
        attribute = ' ' << attribute
      end
      label = label_a + ':' + label_b
      "\n" + indendentation( indent ) + '<' + label + attribute + '>' + value.to_s + '</' + label + '>'
    end

    def indendentation( indent )
      indendentation = String.new
      indent.times do
        indendentation << "\t"
      end
      indendentation
    end

  end # class Hmmer3hit


  class HmmerPerSequenceHit < Hmmer3Hit

    # Sets hit data.
    def initialize(line, line_number)

      # tblout:
      #               tn    tacc      qn    qacc fs_eval fs_scor fs_bias   bst_e bst_scor bst_bias   exp     reg     clu      ov     env      dom     rep     inc   desc
      #                1       2       3       4       5       6       7       8       9      10      11      12      13      14      15       16      17      18     19
      if  line =~ /^(\S*)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(.*)/
        @target_name = $1
        @target_accession = $2
        @query_name = $3
        @query_accession = $4
        @full_sequence_e_value = $5.to_f
        @full_sequence_score = $6.to_f
        @full_sequence_bias =  $7.to_f
        @best_1_domain_e_value = $8.to_f
        @best_1_domain_score = $9.to_f
        @best_1_domain_bias =  $10.to_f
        @domain_number_est_exp = $11.to_i
        @domain_number_est_reg = $12.to_i
        @domain_number_est_clu = $13.to_i
        @domain_number_est_ov  = $14.to_i
        @domain_number_est_env = $15.to_i
        @domain_number_est_dom = $16.to_i
        @domain_number_est_rep = $17.to_i
        @domain_number_est_inc = $18.to_i
        @target_description = $19
      else
        raise ArgumentError, "line "+ line_number.to_s + " is in an unrecognized format [#{line}]"
      end

    end # initialize

    attr_reader :best_1_domain_e_value
    attr_reader :best_1_domain_score
    attr_reader :best_1_domain_bias
    attr_reader :domain_number_est_exp
    attr_reader :domain_number_est_reg
    attr_reader :domain_number_est_clu
    attr_reader :domain_number_est_ov
    attr_reader :domain_number_est_env
    attr_reader :domain_number_est_dom
    attr_reader :domain_number_est_rep
    attr_reader :domain_number_est_inc

    # RDF output

    RDF_HMMER_BEST_1_DOMAIN_E_VALUE = 'best_1_domain_e_value'
    RDF_HMMER_BEST_1_DOMAIN_SCORE = 'best_1_domain_score'
    RDF_HMMER_BEST_1_DOMAIN_BIAS = 'best_1_domain_bias'
    RDF_HMMER_DOMAIN_NUMBER_EST_ESP = 'domain_number_est_exp'
    RDF_HMMER_DOMAIN_NUMBER_EST_REG = 'domain_number_est_reg'
    RDF_HMMER_DOMAIN_NUMBER_EST_CLU = 'domain_number_est_clu'
    RDF_HMMER_DOMAIN_NUMBER_EST_OV = 'domain_number_est_ov'
    RDF_HMMER_DOMAIN_NUMBER_EST_ENV = 'domain_number_est_env'
    RDF_HMMER_DOMAIN_NUMBER_EST_DOM = 'domain_number_est_dom'
    RDF_HMMER_DOMAIN_NUMBER_EST_REP = 'domain_number_est_rep'
    RDF_HMMER_DOMAIN_NUMBER_EST_INC = 'domain_number_est_inc'

    def to_rdf(type = :xml)
      rdf = String.new
      rdf << "\n\t" << '<rdf:Description rdf:about="' + RDF_HMMER_PFAM_RESOURCE + target_accession + '">'
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_TARGET_NAME, target_name, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_QUERY_NAME, query_name, 2)
      if query_accession.length > 1
        rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_QUERY_ACCESSION, query_accession, 2)
      end
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_FULL_SEQUENCE_E_VALUE, full_sequence_e_value, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_FULL_SEQUENCE_SCORE, full_sequence_score, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_FULL_SEQUENCE_BIAS, full_sequence_bias, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_BEST_1_DOMAIN_E_VALUE, best_1_domain_e_value, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_BEST_1_DOMAIN_SCORE, best_1_domain_score, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_BEST_1_DOMAIN_BIAS, best_1_domain_bias, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_NUMBER_EST_ESP, domain_number_est_exp, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_NUMBER_EST_REG, domain_number_est_reg, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_NUMBER_EST_CLU, domain_number_est_clu, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_NUMBER_EST_OV, domain_number_est_ov, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_NUMBER_EST_ENV, domain_number_est_env, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_NUMBER_EST_DOM, domain_number_est_dom, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_NUMBER_EST_REP, domain_number_est_rep, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_NUMBER_EST_INC, domain_number_est_inc, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_TARGET_DESCRIPTION, target_description, 2)
      rdf << "\n" << "\t" << '</rdf:Description> '
      rdf
    end
  end # class HmmerPerSequenceHit

  class HmmerPerDomainHit < Hmmer3Hit

    # Sets hit data.
    def initialize(line, line_number)

      # domtblout:
      #                tn     acc    tlen   query     acc    qlen  Evalue   score    bias      #       of     c-E    i-E     score   bias      hf      ht      af      at     ef      et     acc  desc
      #                 1       2       3       4       5       6       7       8       9      10      11      12     13      14      15       16      17      18      19     20      21      22     23
      if  line =~ /^(\S*)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s*(.*)/
        @target_name = $1
        @target_accession = $2
        @target_length = $3.to_i
        @query_name = $4
        @query_accession = $5
        @query_length = $6.to_i
        @full_sequence_e_value = $7.to_f
        @full_sequence_score = $8.to_f
        @full_sequence_bias =  $9.to_f
        @domain_number = $10.to_i
        @domain_sum = $11.to_i
        @domain_c_e_value = $12.to_f
        @domain_i_e_value = $13.to_f
        @domain_score = $14.to_f
        @domain_bias = $15.to_f
        @hmm_coord_from = $16.to_i
        @hmm_coord_to = $17.to_i
        @ali_coord_from = $18.to_i
        @ali_coord_to = $19.to_i
        @env_coord_from = $20.to_i
        @env_coord_to = $21.to_i
        @acc = $22.to_f
        @target_description = $23
      else
        raise ArgumentError, "line "+ line_number.to_s + " is in a unrecognized format [#{line}]"
      end

    end # initialize

    attr_reader :target_length
    attr_reader :query_length
    attr_reader :domain_number
    attr_reader :domain_sum
    attr_reader :domain_c_e_value
    attr_reader :domain_i_e_value
    attr_reader :domain_score
    attr_reader :domain_bias
    attr_reader :hmm_coord_from
    attr_reader :hmm_coord_to
    attr_reader :ali_coord_from
    attr_reader :ali_coord_to
    attr_reader :env_coord_from
    attr_reader :env_coord_to
    attr_reader :acc

    # RDF output
    RDF_HMMER_TARGET_LENGTH = 'target_length'
    RDF_HMMER_QUERY_LENGTH = 'query_length'
    RDF_HMMER_DOMAIN_NUMBER = 'domain_number'
    RDF_HMMER_DOMAIN_SUM = 'domain_sum'
    RDF_HMMER_DOMAIN_C_E_VALUE = 'domain_c_e_value'
    RDF_HMMER_DOMAIN_I_E_VALUE = 'domain_i_e_value'
    RDF_HMMER_DOMAIN_SCORE = 'domain_score'
    RDF_HMMER_DOMAIN_BIAS = 'domain_bias'
    RDF_HMMER_HMM_COORD_FROM = 'hmm_coord_from'
    RDF_HMMER_HMM_COORD_TO = 'hmm_coord_to'
    RDF_HMMER_ALI_COORD_FROM = 'ali_coord_from'
    RDF_HMMER_ALI_COORD_TO = 'ali_coord_to'
    RDF_HMMER_ENV_COORD_FROM = 'env_coord_from'
    RDF_HMMER_ENV_COORD_TO = 'env_coord_to'
    RDF_HMMER_ACC = 'acc'

    def to_rdf(type = :xml)
      rdf = String.new
      rdf << "\n\t" << '<rdf:Description rdf:about="' + RDF_HMMER_PFAM_RESOURCE + target_accession + '">'
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_TARGET_NAME, target_name, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_TARGET_LENGTH, target_length, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_QUERY_NAME, query_name, 2)
      if query_accession.length > 1
        rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_QUERY_ACCESSION, query_accession, 2)
      end
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_QUERY_LENGTH, query_length, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_FULL_SEQUENCE_E_VALUE, full_sequence_e_value, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_FULL_SEQUENCE_SCORE, full_sequence_score, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_FULL_SEQUENCE_BIAS, full_sequence_bias, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_NUMBER, domain_number, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_SUM, domain_sum, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_C_E_VALUE, domain_c_e_value, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_I_E_VALUE, domain_i_e_value, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_SCORE, domain_score, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_DOMAIN_BIAS, domain_bias, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_HMM_COORD_FROM, hmm_coord_from, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_HMM_COORD_TO, hmm_coord_to, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_ALI_COORD_FROM, ali_coord_from, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_ALI_COORD_TO,ali_coord_to, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_ENV_COORD_FROM, env_coord_from, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_ENV_COORD_TO, env_coord_to, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_ACC, acc, 2)
      rdf << to_rdfxml(RDF_HMMER_HIT, RDF_HMMER_TARGET_DESCRIPTION, target_description, 2)
      rdf << "\n" << "\t" << '</rdf:Description> '
      rdf
    end

  end # class HmmerPerDomainHit

end # module Bio