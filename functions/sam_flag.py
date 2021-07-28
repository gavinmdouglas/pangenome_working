#!/usr/bin/python3

# Taken from https://www.biostars.org/p/13051/
def asbin(n):
    '''converted a number to its binary representative (padded with 0s)'''
    return str(bin(n))[2:].zfill(17)


class human_readable_sam_flags:
    def __init__(self, flag):
        self.read_paired = bool(int(asbin(flag)[-1]))
        self.read_mapped_in_proper_pair = bool(int(asbin(flag)[-2]))
        self.read_unmapped = bool(int(asbin(flag)[-3]))
        self.mate_unmapped = bool(int(asbin(flag)[-4]))
        self.read_reverse_strand = bool(int(asbin(flag)[-5]))
        self.mate_reverse_strand = bool(int(asbin(flag)[-6]))
        self.first_in_pair = bool(int(asbin(flag)[-7]))
        self.second_in_pair = bool(int(asbin(flag)[-8]))
        self.not_primary_alignment = bool(int(asbin(flag)[-9]))
        self.read_fails_platform_vendor_quality_checks = bool(int(asbin(flag)[-10]))
        self.read_is_PCR_or_optical_duplicate = bool(int(asbin(flag)[-11]))
        self.supplementary_alignment = bool(int(asbin(flag)[-12]))


class pertinent_human_readable_sam_flags:
    def __init__(self, flag):
        self.read_paired = bool(int(asbin(flag)[-1]))
        self.read_unmapped = bool(int(asbin(flag)[-3]))
        self.mate_unmapped = bool(int(asbin(flag)[-4]))
        self.first_in_pair = bool(int(asbin(flag)[-7]))
        self.second_in_pair = bool(int(asbin(flag)[-8]))
        self.not_primary_alignment = bool(int(asbin(flag)[-9]))
