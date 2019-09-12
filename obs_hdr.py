"""parse the header lines in an observation batch file"""

import re
from components.obs80.obs80 import parse80

# maps header key to corresponding database table fields
key_col = dict([
    ("COD", "obs_code"), ("CON", "contact_id"), ("NET", "catalog"),
    ("BND", "mag_band"), ("COM", "comments"), ("NUM", "obs_count"),
    ("ACK", "ack_message"), ("AC2", "email"), ("TEL", "tel_id"),
    ("TYP", "comments")
])
ign_col = [("OBS", "batch_people.person_id,role=observer"),
           ("MEA", "batch_people.person_id,role=measurer")]

hdr_re = re.compile(r"({}) (.*)".format('|'.join(key_col.keys())))


def parse(f):
    """reads input file looking for MPC header keys, returns [ {}, ... ]"""
    # submissions can be a series: (head, data)[, ...]
    hs = []
    # head is dict obtained from lines: key SP value NL[, ...]
    # each head starts with key=="COD"
    start_key="COD"
    # data are consecutive 80-character observation records
    def count_obs(obs):
        """count non-header lines that parse as observations"""
        return len([r for r in parse80(obs) if not isinstance(r[0], Exception)])
    obs= []
    hd = None
    for s in f:
        m = hdr_re.match(s) # only matches valid header lines
        if not m:
            if len(s.rstrip())==80: obs.append(s)
            continue
        k, v = m.groups()
        if k == start_key:
            "begin new block"
            if len(hs):
                "close current block"
                hd[key_col['NUM']] = count_obs(obs)
                obs = []
            hd = {}
            hs.append(hd)
        if hd != None:
            c = key_col[k]
            v = v.rstrip()
            hd[c] = v if c not in hd else '\n'.join([hd[c],v])
    # count up last data block
    hd[key_col['NUM']] = count_obs(obs)
    return hs
