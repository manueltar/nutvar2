-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

# patch_66_67_b.sql
#
# Title: Drop stable ID views.
#
# Description:
# The stable ID views, introduced for release 65 as a way of providing a
# degree of backward compatibility, are dropped with this release.

-- Heavens!

DROP VIEW IF EXISTS exon_stable_id;
DROP VIEW IF EXISTS gene_stable_id;
DROP VIEW IF EXISTS operon_stable_id;
DROP VIEW IF EXISTS operon_transcript_stable_id;
DROP VIEW IF EXISTS translation_stable_id;
DROP VIEW IF EXISTS transcript_stable_id;

DROP TABLE IF EXISTS exon_stable_id;
DROP TABLE IF EXISTS gene_stable_id;
DROP TABLE IF EXISTS operon_stable_id;
DROP TABLE IF EXISTS operon_transcript_stable_id;
DROP TABLE IF EXISTS translation_stable_id;
DROP TABLE IF EXISTS transcript_stable_id;



# Patch identifier:
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_66_67_b.sql|drop_stable_id_views');