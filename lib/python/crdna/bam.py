#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

class BamTemplateShim(object):
    def __init__(self, template, keep_comments):
        self.header = template.header
        self.references = template.references
        if not keep_comments and self.header.has_key('CO'):
            self.header = {k:v for k, v in self.header.iteritems() if k != 'CO'}
