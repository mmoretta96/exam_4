"""Tests for the k-mer analyzer starter script."""

from pathlib import Path

import pytest

import kmer_analyzer


def test_validate_sequence_accepts_valid_dna_sequence():
    assert kmer_analyzer.validate_sequence("ATGC", 2) is True


def test_validate_sequence_rejects_sequence_shorter_than_k():
    assert kmer_analyzer.validate_sequence("AT", 3) is False


def test_validate_sequence_rejects_numbers():
    assert kmer_analyzer.validate_sequence("AT1G", 2) is False


def test_validate_sequence_rejects_non_dna_letters():
    assert kmer_analyzer.validate_sequence("ATXG", 2) is False


def test_update_kmer_count_counts_first_observation_once():
    kmer_data = {}

    result = kmer_analyzer.update_kmer_count(kmer_data, "AT", "G")

    assert result["AT"]["count"] == 1
    assert result["AT"]["next_chars"] == {"G": 1}


def test_update_kmer_count_accumulates_repeated_observations():
    kmer_data = {}

    kmer_analyzer.update_kmer_count(kmer_data, "AT", "G")
    result = kmer_analyzer.update_kmer_count(kmer_data, "AT", "G")

    assert result["AT"]["count"] == 2
    assert result["AT"]["next_chars"] == {"G": 2}


def test_count_kmers_with_context_counts_overlapping_kmers():
    result = kmer_analyzer.count_kmers_with_context("ATGA", 2)

    assert result["AT"]["count"] == 1
    assert result["AT"]["next_chars"] == {"G": 1}
    assert result["TG"]["count"] == 1
    assert result["TG"]["next_chars"] == {"A": 1}


def test_count_kmers_with_context_includes_terminal_kmer_total():
    result = kmer_analyzer.count_kmers_with_context("ATGA", 2)

    assert result["GA"]["count"] == 1
    assert result["GA"]["next_chars"] == {}


def test_write_results_to_file_includes_total_counts(tmp_path):
    output_path = tmp_path / "output.txt"
    kmer_data = {
        "AT": {"count": 2, "next_chars": {"G": 2}},
        "TG": {"count": 1, "next_chars": {"A": 1}},
    }

    kmer_analyzer.write_results_to_file(kmer_data, output_path)

    output_text = output_path.read_text()

    assert "AT 2 G:2" in output_text
    assert "TG 1 A:1" in output_text


def test_main_combines_multiple_sequences_before_writing_output(tmp_path, monkeypatch):
    input_path = tmp_path / "fragments.txt"
    output_path = tmp_path / "results.txt"
    input_path.write_text("ATGA\nATGT\n")

    monkeypatch.setattr(
        "sys.argv",
        ["kmer_analyzer.py", str(input_path), "2", str(output_path)],
    )

    kmer_analyzer.main()

    output_text = output_path.read_text()

    assert "AT 2 G:2" in output_text
    assert "TG 2 A:1 T:1" in output_text
