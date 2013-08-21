import unittest
import optparse
import rescale
import pysam
import filecmp

def mock_options(filename,rescale_out,folder):
    """Make the options object with nice values for testing"""
    return optparse.Values({
        "filename":filename,
        "rescale_out":rescale_out,
        "verbose":True,
        "folder":folder,
        "quiet":True
        })


class testPairedFile(unittest.TestCase):
    """Tests if rescaling of a paired end file"""
    def test_paired_end_file(self):
        """Test, paired end SAM file"""
        # here is the example
        # >ref
        # 5p CGA AAA CGA 3p
        # 3p GCT TTT GCT 5p
        # >r001/1
        # CGA
        # >r001/2
        # GCA
        # The sam file looks like this
        #@HD VN:1.5 SO:coordinate
        #@SQ SN:ref LN:9
        #r001 163 ref 1 30 3M = 7 9 CGA III #the normal ones
        #r001 83 ref 7 30 3M = 1 -9 TCG III
        #r002 163 ref 1 30 3M = 7 9 TGA III #With one dam subs
        #r002 83 ref 7 30 3M = 1 -9 CAA III
        #r003 83 ref 1 30 3M = 7 9 TGA III #The reverse complement
        #r003 163 ref 7 30 3M = 1 -9 CAA III
        # 
        #hand calculated the correct rescaled sam file in pe_rescaled_correct.sam
        options = mock_options("rescale_test/pe_test/pe.sam","rescale_test/pe_test/pe_rescaled.sam","rescale_test/pe_test/pe_output/")
        ref = pysam.Fastafile("rescale_test/pe_test/ref.fa")
        rescale.rescale_qual(ref,options,debug=True)
        self.assertTrue(filecmp.cmp("rescale_test/pe_test/pe_rescaled.sam","rescale_test/pe_test/pe_rescaled_correct.sam"))

if  __name__=='__main__':
    unittest.main()
