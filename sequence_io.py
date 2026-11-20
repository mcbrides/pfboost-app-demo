"""Experimental sequence records and deterministic derived FASTA; no PDB lookup."""
import hashlib
import re

AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWYBXZOU')


def checked_sequence(value):
    if not isinstance(value, str) or not value or set(value) - AMINO_ACIDS:
        raise ValueError('Sequence must contain reviewed uppercase amino-acid letters without gaps or whitespace')
    return value


def sequence_records(metadata, entry_id=None):
    """Split only explicitly reviewed separate entities; offsets retain table numbering."""
    identifier = entry_id or metadata.get('dataset_id') or metadata.get('entry_id') or 'sequence'
    if not isinstance(identifier, str) or not re.fullmatch(r'[A-Za-z0-9_.-]+', identifier):
        raise ValueError('Unsafe FASTA record identifier')
    sequence = metadata.get('sequence')
    if sequence in (None, ''):
        return []
    sequence = checked_sequence(sequence)
    review = metadata.get('reviewed_heteromeric_source') or {}
    entities = review.get('entities') or []
    records = []
    if entities:
        normalization = review.get('normalization') or {}
        offsets = normalization.get('entity_offsets_zero_based') or {}
        if normalization.get('global_sequence', sequence) != sequence:
            raise ValueError('Reviewed entity indexing sequence disagrees with metadata.sequence')
        occupied = set()
        for entity in entities:
            name = entity['source_entity_id']
            seq = checked_sequence(entity['sequence'])
            offset = offsets.get(name)
            if type(offset) is not int or offset < 0 or sequence[offset:offset + len(seq)] != seq:
                raise ValueError('Reviewed entity sequence/offset does not match the saved indexing frame')
            positions = set(range(offset, offset + len(seq)))
            if occupied & positions:
                raise ValueError('Overlapping entity indexing ranges')
            occupied |= positions
            suffix = re.sub(r'[^A-Za-z0-9_.-]', '_', name)
            records.append({'id': identifier + '_entity_' + suffix, 'entry_id': identifier,
                            'entity_id': name, 'sequence': seq, 'global_offset_zero_based': offset})
        if occupied != set(range(len(sequence))):
            raise ValueError('Entity records do not cover the saved indexing sequence')
    else:
        records.append({'id': identifier, 'entry_id': identifier,
                        'entity_id': str(metadata.get('entity_id') or ''), 'sequence': sequence,
                        'global_offset_zero_based': 0})
    if len({r['id'] for r in records}) != len(records):
        raise ValueError('FASTA identifiers are not unique')
    for record in records:
        record['sequence_sha256'] = hashlib.sha256(record['sequence'].encode('ascii')).hexdigest()
    return records


def fasta_bytes(records):
    text = ''
    for record in records:
        sequence = checked_sequence(record['sequence'])
        text += '>' + record['id'] + '\n'
        text += '\n'.join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + '\n'
    return text.encode('ascii')
