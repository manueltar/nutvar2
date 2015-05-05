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


# adds a variation set table for structural variation

CREATE TABLE IF NOT EXISTS variation_set_structural_variation (
	structural_variation_id int(10) unsigned NOT NULL,
	variation_set_id int(10) unsigned NOT NULL,
	
	PRIMARY KEY (structural_variation_id,variation_set_id)
);


# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL,'patch', 'patch_64_65_d.sql|adds a variation set table for structural variation');